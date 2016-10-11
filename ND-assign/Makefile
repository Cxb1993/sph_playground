FC=gfortran
FFLAGS=-fopenmp -O3 -Wall -fcheck=all -fdefault-real-8 -fdefault-double-8 #-g -fbacktrace -ffpe-trap=invalid
SRC= kernel.f90 setup.f90 printer.f90 eos.f90 circuit1.f90 circuit2.f90 internal.f90 main.f90
OBJ=$(SRC:.f90=.o)
SUBOBJ=$(addprefix obj/, $(OBJ))

%.o : %.f90
	$(FC) $(FFLAGS) -J mod/ -o obj/$@ -c $<

execute: $(OBJ)
	$(FC) $(FFLAGS) -I mod/ -o $@ $(SUBOBJ)

.PHONY: execute

clean:
	rm -r obj/*
	rm -r mod/*
	rm -r *.png
	rm execute steps/output_*
