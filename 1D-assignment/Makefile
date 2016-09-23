FC=gfortran
FFLAGS=-O3 -Wall -fcheck=all -fdefault-real-8 -fdefault-double-8 #-g -fbacktrace -ffpe-trap=invalid
SRC=main.f90 internal.f90 setup.f90 print.f90 kernel.o
OBJ=$(SRC:.f90=.o)
SUBOBJ=$(addprefix obj/, $(OBJ))

%.o : %.f90
	$(FC) $(FFLAGS) -J mod/ -o obj/$@ -c $<

execute: $(OBJ)
	$(FC) $(FFLAGS) -I mod/ -o $@ $(SUBOBJ)

main.o : internal.o setup.o print.o
internal.o : print.o kernel.o setup.o
.PHONY: execute

clean:
	echo '' > energy.dat
	rm -r obj/*
	rm -r mod/*
	rm execute steps/output_*
