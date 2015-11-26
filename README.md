# mathdws
## A math workspace and libraries for Dyalog APL

This repository contains a math workspace for Dyalog APL with functions for
finding eigenvalues and eigenvectors, matrix inversion and division, and
discrete Fourier transforms.

These functions are implemented as calls into shared libraries. Source code and
build scripts for the libraries are included.

## Building the libraries

The only supported build configuration is cross-compiling from Linux to Windows.

To build (tested on Ubuntu 15.10):

* install some cross-compilers with
`sudo apt-get install gcc-mingw-w64 gfortran-mingw-w64`
* `cd` into the root of the repository
* type `make` (or e.g. `make -j8` if you have 8 CPU cores)

This will generate the release package `math.zip` in the current directory.
