fastppm
--

fastppm is a C++/python libary for fast computation of a maximum likelihood frequency matrix $F \in [0,1]^{m \times n}$ given variant and total read count matrices $$V,D \in \mathbb{N}^{m \times n}$$. This is a key subproblem in cancer phylogeny inference from bulk sequencing data.

## Dependencies

* C++14 compiler
* CMake (>= v3.12)

## Compilation

1. Clone the repository
2. Make sure to initialize all git submodules (`git submodule update --init --recursive`)
3. `mkdir build`
4. `cd build`
5. `cmake ..`
6. `make`

Successful compilation will result in one python library and one executable:

* `fastppm.cpython-310-darwin.so`, which is a Python libary whose extension will vary depending on OS and python libary
* `fastppm-test`, which is an example executable

## Usage

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
