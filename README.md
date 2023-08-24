# SliQEC - A BDD-based Quantum Circuit Equivalence Checker

## Introduction
**SliQEC** is a BDD-based quantum circuit equivalence checker implemented in C/C++ on top of [CUDD](http://web.mit.edu/sage/export/tmp/y/usr/share/doc/polybori/cudd/cuddIntro.html) package.


## Build
Clone the project
```
$ git clone --recurse-submodules https://github.com/joeyweii/sliqec.git
```
Build the checker with CMake, type the command at the root directory.
```
$ mkdir build
$ cd build
$ cmake ..
$ cmake --build .
```

## Execution
The circuit format being checked is `OpenQASM` used by IBM's [Qiskit](https://github.com/Qiskit/qiskit), and our supported gate set now contains Pauli-X, Pauli-Y, Pauli-Z, Hadamard, Phase and its inverse, π/8 and its inverse, Rotation-X with phase π/2, Rotation-Y with phase π/2, Controlled-NOT, Controlled-Z, Toffoli, SWAP, and Fredkin. One can find some example benchmarks in [examples](https://github.com/NTU-ALComLab/SliQEC/tree/main/examples) folder.

The help message concludes the details for execution:

```
$ ./SliQEC --help
Options:
  --help                      produce help message.
  --reorder arg (=1)          allow variable reordering or not.
                              0: disable 1: enable default: 1
  --bitwidth_control arg (=0) bitwidth control when overflowing
                              0: extendBitwidth 1: dropLSB
  --init_bitwidth arg (=4)    initial bitwidth r
                              default: 4
  --circuit1 arg              circuit1 under equivalence checking.

  --circuit2 arg              circuit2 under equivalence checking.

  --approach arg (=miter)     approach of equivalence checking
                              miter: construct miter
                              construct: construct fucntionality
                              simulation: simulation
```
