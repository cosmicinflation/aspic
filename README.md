# ASPIC: Accurate Slow-roll Prediction for Inflationary Cosmology

---

### Summary

This package compiles and install the shared scientific library
[**libaspic**](https://curl.irmp.ucl.ac.be/~chris/aspic.html), a
collection of fast modern fortran routines for computing various
observable quantities used in Cosmology from definite single field
inflationary models.

It aims at providing an efficient, extendable and accurate way of
comparing theoretical inflationary predictions with cosmological
data. As the list of inflationary models is always increasing, you are
encouraged to add support for any model that would not yet be
implemented.

---

### Installation

The code is therefore released as GNU software, you are free and
welcome to expand this code and distribute the source code along the
rules specified in COPYING.

Compilation and installation of the **ASPIC** library

```bash
./configure
  make
  make install
```
See INSTALL for more options.

Checking each models individually (slower)

```bash
  ./configure
  make check
```
In each model subdirectory (located under src/), a test program print
and generate some slow-roll predictions for various values of the reheating
energy density.

Testing the whole library (much slower)

```bash
  ./configure
  make test
```  
This is the equivalent of 'make check' followed by the execution of
all the test programs. The testsuite returns an error if one of them
fails to terminate properly.

---

### Documentation

Checkout the MAN pages for a complete documentation man libaspic

---

### Troubleshooting

Some calls to the "atan()" intrinsic functions use FORTRAN08 support
for complex numbers. This is supported with recent versions of open
source compilers such as gfortran. If you ever encounter an error with
these function calls try to define:
```bash
   export FCFLAGS="-DNOF08"
```  
Then run the standard
```bash
   ./configure"
   make
   make install
```

Some models may require quite extreme fine-tunings according to the
parameter values used. If you get error messages due to numerical
precision limitation, these may be overcome by compiling the library
in quadruple precision (much slower) with

```bash
  ./configure --enable-quad-precision
  make
  make install
```
The new modules are located in $PREFIX/include/aspicq and the library
name is accordingly changed to "libaspicq".

Parallel processing is switched by default using the "-fopenmp"
flag. It is compatible with both gcc, gfortran and other proprietary
compilers. In case of incompatibility with your compiler, parallel
processing can be deactivated with:

```bash
  ./configure --disable-openmp
```

In case you want to specify your own OpenMP compilation flag,
you can do:

```bash
  export FCFLAGS=" -myopenmpflag"
  ./configure --disable-openmp
```

---