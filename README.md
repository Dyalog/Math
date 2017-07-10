# Math

## A math namespace and libraries for Dyalog APL

This repository contains a math namespace for Dyalog APL with functions for
finding eigenvalues, eigenvectors and discrete Fourier transforms.

These functions are implemented as calls into shared libraries. Source code and
build scripts for the libraries are included.

## Using the namespace

Download `math.zip` and extract the files somewhere on your computer. In order
for Windows to find the DLLs you'll have to put them:

* in the current directory (this might be determined by the properties of the
shortcut you use to start Dyalog), or
* in the directory containing the Dyalog executable `dyalog.exe`, or
* somewhere on your `%PATH%`.

## Building the libraries

The only supported build configuration is cross-compiling from Linux to Windows.

To build (tested on Ubuntu 15.10):

* install some cross-compilers with
`sudo apt-get install gcc-mingw-w64 gfortran-mingw-w64`
* `cd` into the root of the repository
* type `make` (or e.g. `make -j8` if you have 8 CPU cores)

This will generate the release package `math.zip` in the current directory.
