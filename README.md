# Aloges library

## Introduction
The **Aloges** library is a collection of Fortran modules providing general-purpose mathematical routines and file /IO utilities used by the *Aloges Lagrangian Model*.<br>
The library has been tested with the **GNU Fortran** (gfortran) **13.3.0** on **Ubuntu 24.04.2 LTS (noble)**.

## Pre-requisistes (Ubuntu/Debian)

### Fortran compiler

Install the compiler toolchain:

  sudo apt update<br>
  sudo apt install gfortran

Check:

gfotran --version

### NetCDF (Scientific data format)
The library includes utilities to handle NetCDF files via **NetCDF-Fortran** bindings: [https://www.unidata.ucar.edu/software/netcdf](https://www.unidata.ucar.edu/software/netcdf)

Install:

sudo apt update<br>
sudo apt install libnetcdff-dev

Optional tools:

sudo apt install netcdf-bin ncview

Verify installation:

nf-config --version<br>
dpkg -l | grep -E "libnetcdf|libnetcdff"

### Plplot (graphic library)
Plotting utilities rely on the **PLplot** Fortran bindings: [https://plplot.sourceforge.net/](https://plplot.sourceforge.net/)

Install (if plotting support is required):

sudo update<br>
sudo apt install libplplot-dev

Recommended device drivers:


sudo apt install plplot-driver-cairo plplot-driver-qt plplot-driver-wxwidgets plplot-driver-xwin plplot-tcl

To disable PLplot support, set the following variable to an empty value in make.inc:


MODULE_PLPLOT=


### LAPACK / BLAS libraries

Several numerical modules require LAPACK and BLAS libraries.

Install:

sudo apt update<br>
sudo apt install gfortran liblapack-dev libblas-dev

Verify installation:

dpkg -l | grep -E "libblas|liblapack"

Typical library installation path on Ubuntu:

/usr/lib/x86_l64-linux-gnu

Example compilation using LAPACK/:

gfortran -O2 -o example example.f90 -llapack -lblas

### Limited-memory BFGS (L-BFGS-B)

The Neural Network module uses the limited-memory BFGS with bounds library, L-BFGS-B. 

Install:

sudo apt update<br>
sudo apt install liblbfgsb-dev

Verify installation:

dpkg - l | grep lbfgsb

Example compilation using L-BFGS-B: 

gfortran -O2 -o example example.f90 -llbfgsb -llapack -lblas

