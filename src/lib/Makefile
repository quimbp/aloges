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

LIBDIR = ../../lib/
INCDIR = ../../include/
INC = $(CDFINC)

TARGET = ../../lib/libaloges.a
OBJ = module_types.o module_constants.o module_dictionary.o         \
      module_nc.o module_math.o module_tools.o module_time.o        \
      module_grid.o module_lineargs.o  module_help.o                \
      module_trajectories.o module_minpack.o module_matrix.o        \
      module_enkf.o module_sysdyn.o module_interp.o                 \
      module_geodetic.o module_codar.o

all: $(OBJ) $(LIBDIR) $(INCDIR) aloges.mod
	rm -f $(TARGET)
	$(INSTALL) *.mod $(INCDIR)
	$(AR) $(TARGET) $(OBJ)

clean:
	rm -f $(TARGET) *.o *.mod 
	if test -d $(INCDIR); then rm -f $(INCDIR)module_types.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_constants.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_dictionary.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_math.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_tools.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_time.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_grid.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_lineargs.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_help.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_trajectories.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_minpack.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_enkf.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_sysdyn.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_interp.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)aloges.mod; fi

$(INCDIR):
	@echo "Directory $(INCDIR) does not exist. Creating it"
	mkdir -p $@

$(LIBDIR):
	@echo "Directory $(LIBDIR) does not exist. Creating it"
	mkdir -p $@

module_types.o: module_types.f90
	$(FC) -c $(FFLAGS) $(INC) module_types.f90

module_constants.o: module_types.o module_constants.f90
	$(FC) -c $(FFLAGS) $(INC) module_constants.f90

module_dictionary.o: module_types.o module_dictionary.f90
	$(FC) -c $(FFLAGS) $(INC) module_dictionary.f90

module_nc.o: module_types.o module_nc.f90
	$(FC) -c $(FFLAGS) $(INC) module_nc.f90

module_math.o: module_types.o module_constants.o module_math.f90
	$(FC) -c $(FFLAGS) $(INC) module_math.f90

module_interp.o: module_types.o module_constants.o module_interp.f90
	$(FC) -c $(FFLAGS) $(INC) module_interp.f90

module_tools.o: module_types.o module_constants.o module_math.o module_tools.f90
	$(FC) -c $(FFLAGS) $(INC) module_tools.f90

module_time.o: module_types.o module_constants.o module_math.o module_tools.o module_time.f90
	$(FC) -c $(FFLAGS) $(INC) module_time.f90

module_matrix.o: module_types.o module_constants.o module_tools.o module_matrix.f90
	$(FC) -c $(FFLAGS) $(INC) module_matrix.f90

module_grid.o: module_types.o module_constants.o module_math.o module_tools.o module_grid.f90
	$(FC) -c $(FFLAGS) $(INC) module_grid.f90

module_lineargs.o: module_types.o module_constants.o module_tools.o module_lineargs.f90
	$(FC) -c $(FFLAGS) $(INC) module_lineargs.f90

module_help.o: module_types.o module_constants.o module_tools.o module_help.f90
	$(FC) -c $(FFLAGS) $(INC) module_help.f90

module_trajectories.o: module_types.o module_constants.o module_tools.o module_trajectories.f90
	$(FC) -c $(FFLAGS) $(INC) module_trajectories.f90

module_minpack.o: module_types.o module_minpack.f90
	$(FC) -c $(FFLAGS) $(INC) module_minpack.f90

module_enkf.o: module_types.o module_constants.o module_math.o module_enkf.f90
	$(FC) -c $(FFLAGS) $(INC) module_enkf.f90

module_sysdyn.o: module_types.o module_constants.o module_math.o module_tools.o module_sysdyn.f90
	$(FC) -c $(FFLAGS) $(INC) module_sysdyn.f90

module_geodetic.o: module_types.o module_tools.o module_math.o module_geodetic.f90
	$(FC) -c $(FFLAGS) $(INC) module_geodetic.f90

module_codar.o: module_types.o module_dictionary.o module_tools.o module_codar.f90
	$(FC) -c $(FFLAGS) $(INC) module_codar.f90

aloges.mod: $(OBJ)
	$(FC) -c $(FFLAGS) $(INC) aloges.f90

