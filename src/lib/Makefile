# ======================================================================== #
# ALOGES PROJECT                                                           #
# Quim Ballabrera, April 2022                                              #
# Institut de Ciencies del Mar, CSIC                                       #
# Created created: 2022-04-14                                              #
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
# ALOGES PROJECT - src/lib Makefile (build static library)                 #
# Last Modified: 2025-10-27                                                #
# -------------------------------------------------------------------------#  

include ../../make.inc
.PHONY: all clean 

LIBDIR := ../../lib/
INCDIR := ../../include/
LIBNAME := libaloges.a
TARGET := $(LIBDIR)/$(LIBNAME)

$(info $(MODULE_PLPLOT))

ifneq ($(strip (MODULE_PLPLOT)),)
  $(info "Compiling with PLPLOT utilities")
  PLSRC = module_plplot.f90
  PLOBJ = module_plplot.o
else
  $(info "Compiling without PLPLOT utilities")
  PLSRC = module_plplot_empty.f90
  PLOBJ = module_plplot_empty.o
endif

OBJ = module_types.o $(PLOBJ) module_constants.o module_color.o     \
      module_tools.o module_timer.o module_dictionary.o module_nc.o \
      module_math.o module_linalg.o module_eofs.o                   \
      module_time.o module_solar.o                                  \
      module_grid.o module_lineargs.o  module_help.o                \
      module_trajectories.o module_minpack.o module_matrix.o        \
      module_enkf.o module_sysdyn.o module_interp.o                 \
      module_geodesy.o module_codar.o module_roms.o                 \
      module_fft.o module_spectra.o module_statistics.o             \
      module_regrid.o

all: $(OBJ) $(LIBDIR) $(INCDIR) aloges.o
	rm -f $(TARGET)
	$(INSTALL) *.mod $(INCDIR)
	$(AR) $(TARGET) $(OBJ) aloges.o

clean:
	rm -f $(TARGET) *.o *.mod 
	if test -d $(INCDIR); then rm -f $(INCDIR)module_types.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_plplot.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_constants.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_color.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_timer.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_dictionary.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_math.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_statistics.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_linalg.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_eofs.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_tools.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_time.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_solar.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_grid.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_lineargs.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_help.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_trajectories.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_minpack.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_enkf.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_sysdyn.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_interp.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_codar.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_geodesy.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_matrix.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_nc.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_roms.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_fft.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_spectra.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)module_regrid.mod; fi
	if test -d $(INCDIR); then rm -f $(INCDIR)aloges.mod; fi

$(INCDIR):
	@echo "Directory $(INCDIR) does not exist. Creating it"
	mkdir -p $@

$(LIBDIR):
	@echo "Directory $(LIBDIR) does not exist. Creating it"
	mkdir -p $@

module_types.o: module_types.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_types.f90

module_plplot.o: module_types.f90 module_plplot.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_plplot.f90

module_plplot_empty.o: module_plplot_empty.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_plplot_empty.f90

module_constants.o: module_types.o module_constants.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_constants.f90

module_color.o: module_color.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_color.f90

module_timer.o: module_types.o module_timer.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_timer.f90

module_dictionary.o: module_types.o module_dictionary.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_dictionary.f90

module_tools.o: module_types.o module_constants.o module_color.o module_tools.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_tools.f90

module_nc.o: module_types.o module_tools.o module_nc.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_nc.f90

module_math.o: module_types.o module_constants.o module_tools.o module_math.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_math.f90

module_regrid.o: module_types.o module_tools.o module_regrid.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_regrid.f90

module_statistics.o: module_types.o module_constants.o module_tools.o module_math.o module_statistics.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_statistics.f90

module_linalg.o: module_types.o module_constants.o module_tools.o module_linalg.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_linalg.f90

module_eofs.o: module_types.o module_constants.o module_tools.o module_linalg.o module_eofs.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_eofs.f90

module_interp.o: module_types.o module_constants.o module_interp.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_interp.f90

module_fft.o: module_types.o module_constants.o module_math.o module_tools.o module_fft.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_fft.f90

module_spectra.o: module_types.o module_constants.o module_math.o module_tools.o module_fft.o module_spectra.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_spectra.f90

module_time.o: module_types.o module_tools.o module_time.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_time.f90

module_solar.o: module_types.o module_tools.o module_time.o module_solar.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_solar.f90

module_matrix.o: module_types.o module_constants.o module_tools.o module_matrix.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_matrix.f90

module_grid.o: module_types.o module_constants.o module_math.o module_tools.o module_grid.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_grid.f90

module_lineargs.o: module_types.o module_constants.o module_tools.o module_lineargs.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_lineargs.f90

module_help.o: module_types.o module_constants.o module_tools.o module_help.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_help.f90

module_trajectories.o: module_types.o module_constants.o module_tools.o module_nc.o module_trajectories.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_trajectories.f90

module_minpack.o: module_types.o module_minpack.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_minpack.f90

module_enkf.o: module_types.o module_constants.o module_math.o module_enkf.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_enkf.f90

module_sysdyn.o: module_types.o module_constants.o module_math.o module_tools.o module_sysdyn.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_sysdyn.f90

module_geodesy.o: module_types.o module_tools.o module_math.o module_geodesy.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_geodesy.f90

module_codar.o: module_types.o module_dictionary.o module_tools.o module_codar.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_codar.f90

module_roms.o: module_types.o module_time.o module_tools.o module_nc.o module_roms.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) module_roms.f90

aloges.o: $(OBJ) aloges.f90
	$(FC) -c $(FFLAGS) $(COMMON_INCS) aloges.f90

