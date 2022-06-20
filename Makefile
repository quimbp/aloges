all:
	(cd src/lib; make)
	(cd src/lagrangian; make)

clean:
	(cd src/lib; make clean)
	(cd src/lagrangian; make clean)
