FC=gfortran
FFLAGS=-fopenmp -O3 -Wall -fcheck=all -fdefault-real-8 -fdefault-double-8 #-g -fbacktrace -ffpe-trap=invalid

SRC=\
	args.f90\
	kernel_cubic.f90\
	kernel_quintic.f90\
	kernel.f90\
	printer.f90\
	eos.f90\
	BC.f90\
	IC.f90\
	circuit1.f90\
	circuit2.f90\
	errcalc.f90\
	bias.f90\
	iterator.f90\
	main.f90
OBJ=$(SRC:.f90=.o)
SUBSRC=$(addprefix src/, $(SRC))
SUBOBJ=$(addprefix obj/, $(OBJ))

%.o : src/%.f90
	@mkdir -p obj
	@mkdir -p output
	$(FC) $(FFLAGS) -J mod/ -o obj/$@ -c $<

execute: $(OBJ)
	$(FC) $(FFLAGS) -I mod/ -o $@ $(SUBOBJ)

hydroshock1:
	make && time ./execute 1 hydroshock 0.1 runresult.info fab
hydroshock2:
	make && time ./execute 2 hydroshock 0.1 runresult.info fab
hydroshock3:
	make && time ./execute 3 hydroshock 0.1 runresult.info fab

infslb1:
	make && time ./execute 1 infslb 0.1 runresult.info fab -1 1 5
infslb2:
	make && time ./execute 2 infslb 0.1 runresult.info fab
infslb3:
	make && time ./execute 3 infslb 0.1 runresult.info fab

hc-sinx1:
	make && time ./execute 3 hc-sinx 0.1 runresult.info fab -1 1 .5

.PHONY: execute

clean:
	rm -rf obj/*; rm -rf mod/*; rm -r *.png; rm execute; rm -f output/*
