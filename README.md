# fastppm: fast perfect phylogeny mixture using tree structured dual dynamic programming

`fastppm` (fast perfect phylogeny mixtures) is a C++/Python libary for 
the fast estimation of the (unknown) frequency matrix $F \in [0,1]^{m \times n}$ 
given variant and total read count matrices $$V,D \in \mathbb{N}^{m \times n}$$
over a fixed n-clonal tree $\mathcal{T}$. `fastppm` provides support for
estimation under a variety of convex loss functions such as the $\ell_1$ loss,
the $\ell_2$ loss, and the negative log binomial loss. `fastppm` can be used
either as a command line tool or as a Python library.

If you find this tool useful in your research, please cite the following paper:

```
```

## Installation

`fastppm` requires depends only on a modern C++14 compiler and CMake. To install, simply clone the 
repository and compile the code using CMake, making sure to initialize all git submodules.

```bash
$ git clone git@github.com:elkebir-group/fastppm.git --recursive
```

To build, run the following commands:

```bash
$ mkdir build
$ cd build
$ cmake ..
$ make
```

The output files consist of the following:
* `fastppm-cli` which is a command line tool.
* `fastppm.cpython-39-x86_64-linux-gnu.so`, which is a Python libary whose extension will vary depending on OS and Python installation.

## Usage

### Option 1: As a command line tool

`fastppm-cli` is a command line tool to regress (estimate) the 
unknown frequency matrix $F$ against the provided variant and total read count 
matrices $V$ and $D$ with respect to a fixed n-clonal tree $\mathcal{T}$. 
The tool requires the following three files as input:
* A file containing the $n$-clonal tree $\mathcal{T}$ in adjacency list format.
* A file containing an $m\text{-by-}n$ variant read count matrix $V$.
* A file containing an $m\text{-by-}n$ total read count matrix $D$.
The nodes in the tree are assumed to be labeled $\\{0, \ldots, n - 1\\}$ and correspond
to the $n$ columns of the variant and total read count matrices. 

The tool has the following usage format:
```
Usage: fastppm [--help] --variant VAR [--version] --total VAR [--weights VAR] --tree VAR --output VAR [--root VAR] [--format VAR] [--loss VAR] [--segments VAR]

Optional arguments:
  -h, --help      shows help message and exits
  -v, --variant   Path to the variant read matrix file [required]
  --version       prints version information and exits
  -d, --total     Path to the total read matrix file [required]
  -w, --weights   Path to the weights matrix file [nargs=0..1] [default: ""]
  -t, --tree      Path to the tree file [required]
  -o, --output    Path to the output file [required]
  -r, --root      Root node of the tree [nargs=0..1] [default: 0]
  -f, --format    Output format, either 'concise' or 'verbose' [nargs=0..1] [default: "concise"]
  -l, --loss      Loss function L_i(.) to use for optimization [nargs=0..1] [default: "l2"]
  -K, --segments  Number of segments, only used when loss function is 'binomial' or 'binomial_K' [nargs=0..1] [default: 10]
 ```

> [!NOTE]
> By default, the tool outputs a JSON file containing only the loss function objective
> and the running time. To include the estimated frequency matrix $F$ and the estimated 
> usage matrix $U$, use the `--format verbose` option.


### Option 2: As a Python library

Ensure the Python library is in the Python path. The following example demonstrates
how to use the Python library to regress the frequency matrix $F$ against the provided
variant and total read count matrices $V$ and $D$ with respect to a fixed n-clonal tree $\mathcal{T}$.
The clonal tree is provided as an adjacency list, and the variant and total read count matrices
are provided as lists of lists. The function `fastppm.regress` returns a dictionary containing
the estimated frequency matrix $F$, the usage matrix $U$, and the loss function objective.

```python
>>> import fastppm
>>> tree  = [[1,2], [3,4], [5,6], [], [], [], []]
>>> var   = [[9, 4, 4, 1, 2, 2, 1]]
>>> tot   = [[10, 10, 10, 10, 10, 10, 10]]
>>> freqs = [[0.9, 0.4, 0.4, 0.1, 0.2, 0.2, 0.1]]
>>> fastppm.regress_counts(tree, var, tot, loss_function="l2")
{'objective': 0.0, 'usage_matrix': [[0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.1]], 'frequency_matrix': [[0.9, 0.4, 0.4, 0.1, 0.2, 0.2, 0.1]]}
>>> fastppm.regress_counts(tree, var, tot, loss_function="binomial")
{'objective': 33.22077103454023, 'usage_matrix': [[0.10001187362880015, 0.10000070598832317, 0.10000070598832315, 0.09999343171700832, 0.19999374434782113, 0.19999374434782113, 0.09999343171700832]], 'frequency_matrix': [[0.8999876377351054, 0.3999878820531526, 0.3999878820531526, 0.09999343171700832, 0.19999374434782113, 0.19999374434782113, 0.09999343171700832]]}
>>> fastppm.regress_frequencies(tree, freqs, loss_function="l2")
{'objective': 0.0, 'usage_matrix': [[0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.1]], 'frequency_matrix': [[0.9, 0.4, 0.4, 0.1, 0.2, 0.2, 0.1]]}
>>> fastppm.regress_counts(tree, var, tot, loss_function="l1")
Traceback (most recent call last):
  File "<stdin>", line 1, in <module>
ValueError: loss_function 'l1' is not yet implemented
```

## Examples

As an example, we first simulate a phylogenetic with $n = 20$ clones and $m = 10$ samples
at a read depth of $c = 10$ with the following command:
```bash
python scripts/simulate.py --mutations 20 --samples 10 --coverage 10 --output examples/sim
```
The files have the following output format:
```bash
$ cat examples/sim_tree.txt
0 5 2 15 17 10
1
2 12
3 14
4
5 13 11
6
7
8 9
9
10
11 8
12 19
13 1
14 6
15 18
16 3
17 4 7
18 16
19
$ cat examples/sim_variant_matrix.txt
10 3 2 1 0 1 3 0 0 0 0 1 4 2 0 5 0 0 1 0
7 0 2 4 0 0 0 0 0 0 0 0 6 0 0 12 2 0 8 2
17 1 0 0 2 4 0 4 3 1 3 3 1 0 0 1 0 2 0 0
1 0 1 2 0 3 2 0 6 0 0 2 0 1 0 4 3 0 2 1
8 3 5 1 0 2 3 0 0 0 0 0 10 1 4 2 2 0 1 3
8 0 1 1 0 7 0 0 2 0 1 3 0 0 0 2 1 1 2 0
15 0 0 1 0 0 0 0 0 0 1 0 0 0 1 2 1 3 0 0
6 0 2 0 0 7 2 0 3 3 0 5 4 0 2 0 1 0 1 1
9 0 1 4 3 0 2 0 0 0 0 0 0 0 1 5 5 6 6 0
9 2 2 1 0 3 0 0 3 0 0 2 4 1 1 10 0 0 3 6
$ cat examples/sim_total_matrix.txt
10 14 10 14 16 8 13 11 7 11 5 18 9 7 5 16 7 6 14 13
7 13 9 10 6 7 6 8 10 8 10 12 13 8 5 12 5 7 13 9
17 9 14 18 15 10 8 18 16 11 9 11 12 7 10 12 6 9 9 6
1 8 15 9 9 7 12 7 9 12 6 8 9 14 7 8 7 4 3 9
8 10 9 10 7 12 9 12 8 11 9 9 12 9 16 8 6 3 9 8
8 11 13 8 14 12 9 8 5 8 9 5 12 10 8 15 12 5 17 12
15 10 7 10 2 14 7 9 12 4 15 10 9 8 8 12 14 6 4 8
6 14 8 5 7 18 6 9 9 11 5 11 12 11 10 9 12 8 5 4
9 11 10 5 8 11 14 7 15 7 9 5 11 16 12 9 10 10 14 4
9 10 13 10 7 10 6 7 14 5 6 10 14 10 12 19 9 8 8 14
```

To run `fastppm` on the simulated data, we then use the following command:
```bash
$ ./build/src/fastppm-cli -v examples/sim_variant_matrix.txt\ 
                          -d examples/sim_total_matrix.txt\ 
                          -t examples/sim_tree.txt\ 
                          -o examples/fastppm_results.json\
                          -f verbose -l l2
```

The output file `examples/fastppm_results.json` contains the resulting objective
value, the estimated frequency matrix $F$, and the estimated usage matrix $U$.
