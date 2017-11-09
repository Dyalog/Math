# `Math` – a namespace of mathematical functions and libraries for Dyalog APL

This repository contains a Math namespace for Dyalog APL with functions for
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

The functions in this namespace include complex arithmetic. Dyalog
represents complex numbers *a+bi* as `aJb`.

### `Eigen`

Monadic function `Eigen` takes an n×n real or complex matrix and returns
an (n+1)×n result of Eigen: `Values⍪⍉↑Vectors`:

    ┌───┬───┬───┬───┐
    │  v a l u e s  │  ── Eigen values
    ├───┼───┼───┼───┤
    │ v │ v │ v │ v │  ┐
    ├ e ┼ e ┼ e ┼ e ┤  │
    │ c │ c │ c │ c │  │
    ├ t ┼ t ┼ t ┼ t ┤  ├─ Eigen vectors.
    │ o │ o │ o │ o │  │
    ├ r ┼ r ┼ r ┼ r ┤  │
    │   │   │   │   │  ┘
    └───┴───┴───┴───┘

`Eigen` has been constructed from LAPACK (Linear Algebra Package) double-
precision C functions which are available as source code from
www.netlib.org/lapack

LAPACK.DLL contains these C functions. They can all be called individually
through `⎕NA`. You need to examining each function's parameters in the
corresponding *.C file (downloaded from the internet) in order to correctly
specify their result and argument types.

For example, look at the `⎕NA` call for `dgeev_` in `Eigen`. Compare this with
the parameters specified in file DGEEV.C . All the other double-precision
real and complex LAPACK functions can be called in this way using `⎕NA`.

Trace the following line in order to see `Eigen` in action:

      test.eigen    ⍝ run and trace 10 times

### `Domino`

`Domino` is equivalent to APL's primitive `⌹` function.

Trace the following line in order to see `Domino` in action:

      test.domino   ⍝ run and trace

### `Fourier`

`Fourier` takes a real or complex array right argument.
The left argument signifies:

 1: Fourier Transform (default).  
¯1: Inverse Fourier Transform.

`Fourier` has been constructed from FFTW (the Fastest Fourier Transform
in the World). The FFTW source code of C functions is available from
www.fftw.org

The FFTW.DLL contains all these functions. They can be called individually
through `⎕NA` after examining each function's parameters in the FFTW
documentation.

Check that

      {⍵=¯1 Fourier Fourier ⍵}↓?((5?5),2)⍴100

Trace the following line in order to see `Fourier` in action:

      test.fourier   ⍝ run and trace

Note that C cover functions dft and idft have been added to the DLL.

These functions are:

```
// dft.c discrete fourier transform

#include <fftw.h>

void dft(int *rank, const int *shape, double *data)
{
   fftwnd_plan plan;
   plan = fftwnd_create_plan(*rank, shape, FFTW_FORWARD, FFTW_IN_PLACE);
   fftwnd_one(plan, (void*)data, 0);
   fftwnd_destroy_plan(plan);
}

// idft.c inverse discrete fourier transform

#include <fftw.h>

void idft(int *rank, const int *shape, double *data)
{
   fftwnd_plan plan;
   plan = fftwnd_create_plan(*rank, shape, FFTW_BACKWARD, FFTW_IN_PLACE);
   fftwnd_one(plan, (void*)data, 0);
   fftwnd_destroy_plan(plan);
}

To see an example of these functions, type:

      test.eigen ⋄ test.domino ⋄ test.fourier
```

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
