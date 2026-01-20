# ======================================================================== #
# ALOGES PROJECT                                                           #
# Quim Ballabrera, April 2022                                              #
# Institut de Ciencies del Mar, CSIC                                       #
# Last Modified: 2022-04-14                                                #
#                                                                          #
# Copyright (C) 2022, Joaquim Ballabrera                                   #
#                                                                          #
# This program is free software: you can redistribute it and/or modify     #
# it under the terms of the GNU Lesser General Public License as published #
# by the Free Software Foundation, either version 3 of the License, or     #
# (at your option) any later version.                                      #
#                                                                          #
# This program is distributed in the hope that it will be useful,          #
# but WITHOUT ANY WARRANTY; without even the implied warranty of           #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.                     #
# See the Lesser GNU General Public License for more details.              #
#                                                                          #
# You should have received a copy of the GNU Lesser General                #
# Public License along with this program.                                  #
# If not, see <http://www.gnu.org/licenses/>.                              #
# -------------------------------------------------------------------------#
# Build configuration (portable)                                           #
# Compiler options and libraries
# Last Modified: 2025-10-27                                                #
# -------------------------------------------------------------------------#  

# Modify according operative system and library paths:
#
# OS          = Ubuntu or Centos
# FC          = Fortran95 compiler
# BUILD       = release or debug
# WITH_PLPLOT = Switch to use PLPLOT or not (1-use, 0-don't use)
# FFLAGS      = Compiler options
# INSTALL     = Program used to copy files setting theis mode bits, etc.
# AR          = Library archive manager 
# NCDF_INC    = Path to the NetCDF Fortran 90 modules
# NCDF_LIB    = Path to the NetCDF Fortran 90 libraries
# PLPLOT_INC  = Path to the PLPLOT Fortran 90 modules
# PLPLOT_LIB  = Path to the PLPLOT Fortran 90 libraries
# LBFGSB_LIB  = Library links to LBFGSB/LAPACK/BLAS
# 

# OS and compiler options:
OS = Ubuntu
FC = gfortran

# Build type: release or debug  
BUILD = debug
BUILD = release

# Feature toggles
WITH_PLPLOT = 1

# Installation tool and archiver  
# AR = ar cq
INSTALL = install  
AR      = ar rcs 

# -------------------------------------------------------------------------#  
# Compilation flags
# -------------------------------------------------------------------------#  
  
# Build-independent compilation flags
BASE_FFLAGS := -ffree-line-length-none  

# Debug vs Release flags  
ifeq ($(BUILD),debug)  
  FFLAGS = $(BASE_FFLAGS) -g -O0 -fcheck=all -fbounds-check \
     -fbacktrace -Wall -Wextra -Wno-maybe-uninitialized     \
     -ffpe-summary=zero,invalid,overflow,underflow  
else  
# Keep defaults portable; allow opting into native/unroll via EXTRA_OPT  
  FFLAGS = $(BASE_FFLAGS) -O3  
# Example opt-ins (uncomment or pass via env if desired):  
# EXTRA_OPT ?= -march=native -funroll-loops -fomit-frame-pointer  
  EXTRA_OPT = 
  FFLAGS += $(EXTRA_OPT)  
endif  
  

# -------------------------------------------------------------------------#  
# NetCDF             
# -------------------------------------------------------------------------#  
# dpkg -L libnetcdff-dev
  
ifndef NCDF_INC
  NCDF_INC = -I$(shell nf-config --includedir)
endif  

ifndef NCDF_LIB
  NCDF_LIB = $(shell nf-config --flibs)
endif

# -------------------------------------------------------------------------#  
# PLPLOT (optional)  
# -------------------------------------------------------------------------#  
# dpkg -L libplplot-dev
  
ifeq ($(WITH_PLPLOT),1)  
  # Allow overrides; provide Ubuntu-friendly defaults  
  PLPLOT_INC ?= -I/usr/lib/x86_64-linux-gnu/fortran/modules/plplot/  
  PLPLOT_LIB ?= -lplplotfortran  
  MODULE_PLPLOT ?= module_plplot.o  
else  
  PLPLOT_INC :=  
  PLPLOT_LIB :=  
  MODULE_PLPLOT :=  
endif  

# -------------------------------------------------------------------------#  
# LBFGSB/LAPACK/BLAS
# -------------------------------------------------------------------------#  
# dpkg -L liblbfgsb-dev liblapack-dev libblas-dev
LBFGSB_LIB = -llbfgsb -llapack -lblas


# -------------------------------------------------------------------------#  
# Adding all up together
# -------------------------------------------------------------------------#  
# Common compile-time includes (used by subdir makefiles)  

COMMON_INCS := $(NCDF_INC) $(PLPLOT_INC)  
  
# Link-time libs (used by subdir makefiles)  
COMMON_LIBS := $(NCDF_LIB) $(PLPLOT_LIB)  $(LBFGSB_LIB)
  
# Default install prefix  
# The variable ALOGES, defined in .bashrc
#PREFIX ?= /usr/local  
PREFIX = $(ALOGES)
