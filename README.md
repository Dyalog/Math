# Math

## A math namespace and libraries for Dyalog APL

This repository contains a math namespace for Dyalog APL with functions for
finding eigenvalues, eigenvectors and discrete Fourier transforms.

These functions are implemented as calls into shared libraries. Source code and
build scripts for the libraries are included.

## Using the namespace

Download an appropriate `Math-*.zip` from the releases page and extract the
files somewhere on your computer.

In order for Windows to find the DLLs you'll have to put them:

* in the current directory (this might be determined by the properties of the
shortcut you use to start Dyalog), or
* in the directory containing the Dyalog executable `dyalog.exe`, or
* somewhere on your `%PATH%`.

## Building the libraries

The supported build configurations are:
* Linux (x86 and x86-64)
* Cross-compiling from Linux to Windows (x86 and x86-64)

To build natively:

* install some 32- and 64-bit Fortran compilers with
`sudo apt-get install gfortran-multilib`
* `cd` into the root of the repository
* type `make -f Makefile.linux` to make `Math-linux.zip`

To cross-compile (tested on Ubuntu 17.04):

* install some cross-compilers with
`sudo apt-get install gfortran-mingw-w64`
* `cd` into the root of the repository
* type `make -f Makefile.windows` to make `Math-windows.zip`
