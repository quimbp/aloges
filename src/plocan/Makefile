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

TARGETS = plocan_merge

all: $(BINDIR) $(TARGETS)

$(BINDIR):
	@echo "Directory $(BINDIR) does not exist. Creating it"
	mkdir -p $@

it:
	make clean
	make all

clean:
	rm -f $(TARGETS) 
	(cd $(BINDIR); rm -f $(TARGETS))

plocan_merge: plocan_merge.f90
	$(FC) $(FFLAGS) -o plocan_merge plocan_merge.f90 $(INC) $(LIB)
	$(INSTALL) plocan_merge $(BINDIR)

