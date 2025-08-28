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
#
# Compiler options and libraries
# The code has been created and tested using gfortran.
# The NF90_LIB and NF90_INC variables had been fed from the
# bash environment.

# Modify according operative system and library paths:
#
# OS       = Ubuntu or Centos
# FC       = Fortran95 compiler
# FFLAGS   = Compiler options
# INSTALL  = Program used to copy files setting theis mode bits, etc.
# AR       = Library archive manager 
# NF90_INC = Path to the  NetCDF Fortran 90 modules
# NF90_LIB = Path to the  NetCDF Fortran 90 libraries
# PLPLOT_INC = Path to the PLPLOT Fortran 90 modules
# PLPLOT_LIB = Path to the PLPLOT Fortran 90 libraries
# 
OS = Ubuntu
FC = gfortran
FFLAGS = -O3 
FFLAGS = -fbounds-check
FFLAGS = -fbounds-check -ffpe-summary=zero,invalid,overflow,underflow -ffree-line-length-none
FFLAGS = -ffpe-summary=zero,invalid,overflow,underflow -ffree-line-length-none

INSTALL = install
AR = ar cq

# Other configurations might need manual specification
#
ifeq ($(OS), Centos)
  # Centos:
  NF90_INC = /usr/lib64/gfortran/modules
  NF90_LIB = /usr/lib64
else
  # Ubuntu:
  NF90_INC = $(shell nc-config --includedir)
  NF90_LIB = $(shell nc-config --libdir)
endif

CDFINC = -I$(NF90_INC)
CDFLIB = -L$(NF90_LIB) -lnetcdff

# PLPLOT paths for ubuntu
#
PLPLOT_INC = -I/usr/lib/x86_64-linux-gnu/fortran/modules/plplot/
PLPLOT_LIB = -lplplotfortran
