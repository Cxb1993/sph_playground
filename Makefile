FC=gfortran
FFLAGS=-O3 -Wall -fcheck=all -fdefault-real-8 -fdefault-double-8 #-g -fbacktrace -ffpe-trap=invalid
SRC=main.f90 internal.f90 setup.f90 print.f90
OBJ=$(SRC:.f90=.o)

%.o : %.f90
	$(FC) $(FFLAGS) -o $@ -c $<

execute: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ)

main.o : internal.o setup.o print.o

.PHONY: execute

clean:
	rm *.o *.mod *.dat execute fort.18 step_*
