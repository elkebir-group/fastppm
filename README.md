# fastppm: fast perfect phylogeny mixture using tree structured dual dynamic programming

`fastppm` (fast perfect phylogeny mixtures) is a C++/Python libary for 
the fast estimation of the frequency matrix $F \in [0,1]^{m \times n}$ 
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
* `fastppm.cpython-310-darwin.so`, which is a Python libary whose extension will vary depending on OS and python libary

## Usage

### Option 1: As a command line tool

`fastppm-cli` is a command line tool that can be used to estimate the 
frequency matrix $F$ given the variant and total read count 
matrices $V$ and $D$ over a fixed n-clonal tree $\mathcal{T}$. 
The tool requires the following three inputs:
* A file containing the $n$-clonal tree $\mathcal{T}$ in adjacency list format.
* A file containing an $m$-by-$n$ variant read count matrix $V$.
* A file containing an $m$-by-$n$ total read count matrix $D$.
The nodes in the tree are assumed to be labeled from $0$ to $n-1$ and correspond
to the columns of the variant and total read count matrices. 

The tool outputs a JSON file containing the estimated frequency matrix $F$,
the usage matrix $U$, the loss function objective, and the running time.

The tool has the following usage format:
```
Usage: fastppm [--help] --variant VAR [--version] --total VAR --tree VAR --output VAR [--root VAR] [--loss VAR]

Optional arguments:
  -h, --help     shows help message and exits
  -v, --variant  Path to the variant read matrix file [required]
  --version      prints version information and exits
  -d, --total    Path to the total read matrix file [required]
  -t, --tree     Path to the tree file [required]
  -o, --output   Path to the output file [required]
  -r, --root     Root node of the tree [nargs=0..1] [default: 0]
  -l, --loss     Loss function L_i(.) to use for optimization [nargs=0..1] [default: "binomial"]
 ```

### Option 2: As a Python library

```
import fastppm

a = fastppm.CVXTree(7)

tree = [[1,2],[3,4],[5,6],[],[],[],[]]
root = 0

var = [9,4,4,1,2,2,1]
ref = [1,6,6,9,8,8,8]

ll = a.optimize(tree, root, var, ref)

print("Negative log likelihood:", ll)
print("Frequency matrix:", a.F())
```

## Examples

As an example, we first simulate a phylogenetic with 10 clones and 3 samples
at a read depth of 10 with the following command:
```bash
$ python scripts/simulation.py --mutations 10 --clones 10 --samples 3 --coverage 10 --seed 6 --output examples/sim
```
The files have the followign output format:
Then, we estimate the frequency matrix using the following command:
```bash
$ cat examples/sim_tree.txt
0 1 8 4
1 5
2
3
4
5
6
7 9
8 2 7 6
9 3
$ cat examples/sim_collapsed_variant_matrix.txt
9 9 0 3 0 8 0 3 6 6
11 5 0 0 5 0 0 1 1 0
13 3 0 0 7 0 0 0 0 0
$ cat examples/sim_collapsed_total_matrix.txt
9 14 7 9 11 12 10 10 10 8
11 8 9 10 12 10 9 7 11 9
13 8 8 13 15 9 15 11 6 7
```

To run `fastppm` on the simulated data, we use the following command:
```bash
$ ./build/src/fastppm-cli -v examples/sim_collapsed_variant_matrix.txt\ 
                          -d examples/sim_collapsed_total_matrix.txt\ 
                          -t examples/sim_tree.txt\ 
                          -o examples/sim_results.json
```

The output file `examples/sim_results.json` contains the output,
which should match the simulated frequency matrix `examples/sim_obs_frequency_matrix.txt`.
