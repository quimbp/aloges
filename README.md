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

gfortran --version

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

By default, the library will be compiled with plotting capabilities. In the file **make.inc** file, the following statement can be found:

WITH_PLPLOT = 1

To disable PLplot support, set in **make.inc**:

WITH_PLPLOT = 0


### LAPACK / BLAS libraries

Several numerical modules require LAPACK and BLAS libraries.

Install:

sudo apt update<br>
sudo apt install liblapack-dev libblas-dev

Verify installation:

dpkg -l | grep -E "libblas|liblapack"

Typical library installation path on Ubuntu:

/usr/lib/x86_l64-linux-gnu

Example compilation using LAPACK/BLAS:

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

## Install

The installation of the **Aloges** library and applications has three steps: Download the code, compile and system setup:

### Download the code

To download the code from the **git** repository:

git clone https://github.com/quimbp/aloges.git

### Compile the code

cd aloges

*Edit, if necessary, the file **make.inc**. In this file the user can change the compiler, the compiler options, library paths and linking options. In this file, the variable **WITH_PLPLOT** can be used to activate (default) and deactivate the library's plotting capabilities.*

Compilation:

make

A *script* is available to compile a Fortran program using the options defined in the file **make.inc**. To generate such a script:

make afort

*It will create a bash script in the *bin/* folder that can be used to compile your own applications using the **Aloges** library:

afort your_program.f90


### Setting up the linux environment

Inside the **.bashrc** file in your home folder, define the following variable:

export ALOGES=SOMEPATH/aloges

Then, add:

PATH="$ALOGES/bin:$PATH"

In this way, your sistem will be able to locate the applications and scripts present in your folder *SOMEPATH/aloges/bin/*

## Examples

A series of examples using various modules of the **Aloges** library can be found in the folder src/examples

cd src/examples

To compile all the examples:

make

Each example can be individually compiled using the script **afort**:

afort example_constants.f90

or 

afort example_constants
