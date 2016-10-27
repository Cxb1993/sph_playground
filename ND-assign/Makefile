FC=gfortran
FFLAGS=-fopenmp -O3 -Wall -fcheck=all -fdefault-real-8 -fdefault-double-8 -fstack-arrays #-g -fbacktrace -ffpe-trap=invalid
SRC=kernel.f90\
		printer.f90\
		eos.f90\
		internal/BC.f90\
		internal/IC.f90\
		internal/circuit1.f90\
		internal/circuit2.f90\
		internal/errcalc.f90\
		internal.f90\
		main.f90
OBJ=$(SRC:.f90=.o)
SUBOBJ=$(addprefix obj/, $(OBJ))

%.o : %.f90
	@mkdir -p obj/internal
	@mkdir -p obj/tempr
	@mkdir -p output
	$(FC) $(FFLAGS) -J mod/ -o obj/$@ -c $<

execute: $(OBJ)
	$(FC) $(FFLAGS) -I mod/ -o $@ $(SUBOBJ)

hydroshock1:
	make && time ./execute 1 hydroshock
hydroshock2:
	make && time ./execute 2 hydroshock
hydroshock3:
	make && time ./execute 3 hydroshock

tempr1:
	make && time ./execute 1 heatslab
tempr2:
	make && time ./execute 2 heatslab
tempr3:
	make && time ./execute 3 heatslab

.PHONY: execute

clean:
	rm -rf obj/*; rm -rf mod/*; rm -r *.png; rm execute; rm -f output/*
