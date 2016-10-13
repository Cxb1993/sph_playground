FC=gfortran
FFLAGS=-fopenmp -O3 -Wall -fcheck=all -fdefault-real-8 -fdefault-double-8 #-g -fbacktrace -ffpe-trap=invalid
SRC=kernel.f90 printer.f90 eos.f90 purehydro/setup.f90 purehydro/circuit1.f90 purehydro/circuit2.f90 tempr/setup.f90 tempr/circuit1.f90 tempr/circuit2.f90 internal.f90 main.f90
OBJ=$(SRC:.f90=.o)
SUBOBJ=$(addprefix obj/, $(OBJ))

%.o : %.f90
	@mkdir -p obj/purehydro
	@mkdir -p obj/tempr
	$(FC) $(FFLAGS) -J mod/ -o obj/$@ -c $<

execute: $(OBJ)
	$(FC) $(FFLAGS) -I mod/ -o $@ $(SUBOBJ)

hydroshock1:
	make && time ./execute 1 purehydroshock
hydroshock2:
	make && time ./execute 2 purehydroshock
hydroshock3:
	make && time ./execute 3 purehydroshock

tempr1:
	make && time ./execute 1 temperhomog01

.PHONY: execute

clean:
	rm -rf obj/*; rm -rf mod/*; rm -r *.png; rm execute; rm -f steps/output_*
