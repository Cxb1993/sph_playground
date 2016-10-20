FC=gfortran
FFLAGS=-fopenmp -O3 -Wall -fcheck=all -fdefault-real-8 -fdefault-double-8 #-g -fbacktrace -ffpe-trap=invalid
SRC=kernel.f90\
		printer.f90\
		eos.f90\
		internal/BC.f90\
		internal/circuit1.f90\
		internal/circuit2.f90\
		internal/IC_hydro_shock.f90\
		internal/IC_heat_cond.f90\
		internal.f90\
		main.f90
OBJ=$(SRC:.f90=.o)
SUBOBJ=$(addprefix obj/, $(OBJ))

%.o : %.f90
	@mkdir -p obj/internal
	@mkdir -p obj/tempr
	@mkdir -p steps
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
	make && time ./execute 1 temperhomog01
tempr2:
	make && time ./execute 2 temperhomog01

.PHONY: execute

clean:
	rm -rf obj/*; rm -rf mod/*; rm -r *.png; rm execute; rm -f steps/*
