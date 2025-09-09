# aloges model and library
Aloges Lagrangian Model

## Introduction
The aloges model and library contains a set of Fortran fuctions solving some general-purpose matrhomatical and file processing used by the Aloges Lagrangian Model.
This set of functions have been tested with the gfortran compiler (13.3.0) within an Ubuntu system (Ubuntu 24.04.2 LTS, Release: 24.04, noble).


## Pre-requisistes

The aloges library has been compiled and tested with the gfortran compiler. In Ubuntu, to install the compiler:

> sudo apt-get update

> sudo apt-get install gfortran


The library includes plotting utilities built upon the PLPLOT library (https://plplot.sourceforge.net/) bindings for the Fortran language.
Before compiling the code, if these utilities are wanted, the user must install the fortran library:

> sudo apt-get install libplplotfortran0
>
> sudo apt-get install libplplot-dev
>
> sudo apt-get install plplot-driver-cairo plplot-driver-qt plplot-driver-wxwidgets plplot-driver-xwin plplot-tcl

To skip the plotting utilities of the aloges library, leave empty the variable

> MODULE_PLPLOT=

in the file make.inc

