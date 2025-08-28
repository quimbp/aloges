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
# To compile the library and lagrangian model:
# > make
# To configure the script afort, use:
# > make afort

include make.inc

all:
	(cd src/lib; make)
	(cd src/lagrangian; make)

clean:
	(cd src/lib; make clean)
	(cd src/lagrangian; make clean)
	(rm -f bin/afort)

afort:
	(cd scripts; sed 's@#ALOGES_ROOT.*@ALOGES_ROOT=${PWD}@' afort.template > afort.tmp1)
	(cd scripts; sed 's@#FC.*@FC=${FC}@' afort.tmp1 > afort.tmp2)
	(cd scripts; sed 's@#FFLAGS.*@FFLAGS="'"${FFLAGS}"'"@' afort.tmp2 > afort.tmp3)
	(cd scripts; sed 's@#CDFINC.*@CDFINC="'"${CDFINC}"'"@' afort.tmp3 > afort.tmp4)
	(cd scripts; sed 's@#CDFLIB.*@CDFLIB="'"${CDFLIB}"'"@' afort.tmp4 > afort.tmp5)
	(cd scripts; sed 's@#PLPLOT_INC.*@PLPLOT_INC="'"${PLPLOT_INC}"'"@' afort.tmp5 > afort.tmp6)
	(cd scripts; sed 's@#PLPLOT_LIB.*@PLPLOT_LIB="'"${PLPLOT_LIB}"'"@' afort.tmp6 > afort.tmp7)
	(mv scripts/afort.tmp7 bin/afort; chmod +x bin/afort; rm -f scripts/afort.tmp*)
