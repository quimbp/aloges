include make.inc

all:
	(cd src/lib; make)
	(cd src/lagrangian; make)
	make afort

clean:
	(cd src/lib; make clean)
	(cd src/lagrangian; make clean)
	(rm -f bin/afort)

afort:
	#(cd scripts;sed 's@#ALOGES_ROOT.*@ALOGES_ROOT="'"${PWD}"'"@' afort.template > ../bin/afort; chmod +x ../bin/afort)
	(cd scripts; sed 's@#ALOGES_ROOT.*@ALOGES_ROOT=${PWD}@' afort.template > afort.tmp1)
	(cd scripts; sed 's@#FC.*@FC=${FC}@' afort.tmp1 > afort.tmp2)
	(cd scripts; sed 's@#FFLAGS.*@FFLAGS="'"${FFLAGS}"'"@' afort.tmp2 > afort.tmp3)
	(cd scripts; sed 's@#CDFINC.*@CDFINC="'"${CDFINC}"'"@' afort.tmp3 > afort.tmp4)
	(cd scripts; sed 's@#CDFLIB.*@CDFLIB="'"${CDFLIB}"'"@' afort.tmp4 > afort.tmp5)
	(mv scripts/afort.tmp5 bin/afort; chmod +x bin/afort; rm -f scripts/afort.tmp*)
