# Aloges model and library
Aloges Lagrangian Model

## Introduction
The aloges model and library contains a set of Fortran fuctions solving some general-purpose matrhomatical and file processing used by the Aloges Lagrangian Model.
This set of functions have been tested with the gfortran compiler (13.3.0) within an Ubuntu system (Ubuntu 24.04.2 LTS, Release: 24.04, noble).

## Pre-requisistes

### Fortran compiler

The aloges library has been compiled and tested with the gfortran compiler. In Ubuntu, to install the compiler:

> sudo apt-get update

> sudo apt-get install gfortran

### NetCDF (Scientific data format)
The library includes utilities to handle NetCDF files. These are buil upon the NetCDF (https://www.unidata.ucar.edu/software/netcdf) bindings for the Fortran language. Before compiling the code the user must install them:

> sudo apt-get install libnetcdff-dev

It is also advisable to install the following packages:

> sudo apt-get install netcdf-bin ncview

### Plplot (graphic library)
The library includes plotting utilities built upon the PLPLOT library (https://plplot.sourceforge.net/) bindings for the Fortran language.
Before compiling the code, if these utilities are wanted, the user must install the fortran library:

> sudo apt-get install libplplot-dev

The following device drivers can be of help:

> sudo apt-get install plplot-driver-cairo plplot-driver-qt plplot-driver-wxwidgets plplot-driver-xwin plplot-tcl

To skip the plotting utilities of the aloges library, leave empty the variable

> MODULE_PLPLOT=

in the file make.inc

### Lapack / Blas libraries

In annticipation of using Lapack and BLAS libraries to accelerate the code,
the following libraries should be also installed.

> sudo apt update
> sudo apt install gfortran libblas-dev liblapack-dev

