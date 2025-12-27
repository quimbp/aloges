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
# OS       = Ubuntu or Centos
# FC       = Fortran95 compiler
# BUILD    = release or debug
# PLPLOT   = Switch to use PLPLOT or not (1-use, 0-don't use)
# FFLAGS   = Compiler options
# INSTALL  = Program used to copy files setting theis mode bits, etc.
# AR       = Library archive manager 
# NF90_INC = Path to the  NetCDF Fortran 90 modules
# NF90_LIB = Path to the  NetCDF Fortran 90 libraries
# PLPLOT_INC = Path to the PLPLOT Fortran 90 modules
# PLPLOT_LIB = Path to the PLPLOT Fortran 90 libraries
# LAPACK_LIB = Path to the Lapack Fortran 90 libraries
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
  
# Other configurations might need manual specification
#
#ifeq ($(OS), Centos)
#  # Centos:
#  NF90_INC = /usr/lib64/gfortran/modules
#  NF90_LIB = /usr/lib64
#else
#  # Ubuntu:
#  NF90_INC = $(shell nc-config --includedir)
#  NF90_LIB = $(shell nc-config --libdir)
#endif
#CDF_INC = -I$(NF90_INC)
#CDF_LIB = -L$(NF90_LIB) -lnetcdff
# Allow users to predefine CDF_INC/CDF_LIB (highest priority)  
# Otherwise try nf-config (Fortran-specific)  

ifndef CDF_INC  
  NF_CONFIG := $(shell command -v nf-config 2>/dev/null)  
  ifneq ($(NF_CONFIG),)  
    NF90_INC := $(shell nf-config --includedir 2>/dev/null)  
    CDF_INC   := -I$(NF90_INC)  
  else  
    # Try nc-config (C + Fortran libs discovery)  
    NC_CONFIG := $(shell command -v nc-config 2>/dev/null)  
    ifneq ($(NC_CONFIG),)  
      NF90_INC := $(shell nc-config --includedir 2>/dev/null)  
      CDF_INC   := -I$(NF90_INC)  
    else  
      # Fallback (user must ensure paths are visible to the compiler)  
      CDF_INC   ?=  
    endif  
  endif  
endif  

ifndef CDF_LIB
  ifneq ($(NF_CONFIG),)
    # nf-config --flibs returns full link flags (recommended)
    CDF_LIB := $(shell nf-config --flibs 2>/dev/null)
  else ifneq ($(NC_CONFIG),)
    NF90_LIBDIR := $(shell nc-config --libdir 2>/dev/null)
    # Conservative default if only nc-config is available
    CDF_LIB := -L$(NF90_LIBDIR) -lnetcdff
  else
    # Fallback: rely on system linker paths
    CDF_LIB ?= -lnetcdff
  endif
endif

# -------------------------------------------------------------------------#  
# PLPLOT (optional)  
# -------------------------------------------------------------------------#  
  
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

## PLPLOT paths for ubuntu
## To graphic module module_plplot requires the PLPLOT package
## In ubuntu systems:
## sudo apt install plplot-doc plplot-dev 
## To avoid compiling this module, let it empty
##MODULE_PLPLOT = 
#MODULE_PLPLOT = module_plplot.o
#PLPLOT_INC = -I/usr/lib/x86_64-linux-gnu/fortran/modules/plplot/
#PLPLOT_LIB = -lplplotfortran

# -------------------------------------------------------------------------#  
# LAPACK (optional)  
# -------------------------------------------------------------------------#  
LAPACK_LIB = -llapack -lblas
  
# Common compile-time includes (used by subdir makefiles)  
COMMON_INCS := $(CDF_INC) $(PLPLOT_INC)  
  
# Link-time libs (used by subdir makefiles)  
COMMON_LIBS := $(CDF_LIB) $(PLPLOT_LIB)  $(LAPACK_LIB)
  
# Default install prefix  
# The variable ALOGES, defined in .bashrc
#PREFIX ?= /usr/local  
PREFIX = $(ALOGES)
