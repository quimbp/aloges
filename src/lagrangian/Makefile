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

include ../../make.inc

#DBFLAGS = -ffpe-summary=zero,invalid,overflow,underflow
#DBFLAGS = -ffpe-trap=zero,invalid,overflow,underflow
DBFLAGS =
FCFLAGS = $(FFLAGS) $(DBFLAGS)

BINDIR = ../../bin/
INCDIR = ../../include/
LIBDIR = ../../lib/

INC = -I$(INCDIR) $(CDFINC)
LIB = -L$(LIBDIR) -laloges $(CDFLIB)

TARGET = alm
OBJ = module_alm.o module_forcing.o module_float.o module_model.o module_options.o \
      module_fitting.o

all: $(BINDIR) $(TARGET) $(OBJ) main.f90

$(BINDIR):
	@echo "Directory $(BINDIR) does not exist. Creating it"
	mkdir -p $@

it:
	make clean
	make all

clean:
	rm -f $(TARGET) $(OBJ) *.mod 
	if test -d $(BINDIR); then rm -f $(BINDIR)alm; fi

$(TARGET): $(OBJ) main.f90
	rm -f $(TARGET)
	$(FC) $(FFLAGS) -o $(TARGET) main.f90 $(OBJ) $(INC) $(LIB)
	$(INSTALL) $(TARGET) $(BINDIR)

module_alm.o: module_alm.f90
	$(FC) -c $(FFLAGS) $(INC) module_alm.f90

module_forcing.o: module_forcing.f90
	$(FC) -c $(FFLAGS) $(INC) module_forcing.f90

module_float.o: module_float.f90
	$(FC) -c $(FFLAGS) $(INC) module_float.f90

module_model.o: module_model.f90
	$(FC) -c $(FFLAGS) $(INC) module_model.f90

module_fitting.o: module_fitting.f90
	$(FC) -c $(FFLAGS) $(INC) module_fitting.f90

module_options.o: module_forcing.o module_float.o module_model.o module_fitting.o module_options.f90
	$(FC) -c $(FFLAGS) $(INC) module_options.f90
