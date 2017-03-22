FC=gfortran
FFLAGS=-fdefault-real-8 -fdefault-double-8\
 				-Ofast\
				-g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace\
				-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow,invalid -finit-real=nan\
				-fopenmp

SRC=\
	const.f90\
	kernel/base_cubic.f90\
	kernel/base_quintic.f90\
	kernel/base_external.f90\
	kernel/base_gaus.f90\
	kernel/base_sinq.f90\
	kernel/n2movedgaus.f90\
	kernel/n2ext.f90\
	kernel/kernel.f90\
	args.f90\
	printer.f90\
	eos.f90\
	BC.f90\
	utils.f90\
	IC/initpositions.f90\
	IC/semiuniform.f90\
	IC/IC.f90\
	neighboursearch.f90\
	circuit1.f90\
	circuit2.f90\
	errcalc.f90\
	errteylor.f90\
	iterator.f90\
	main.f90

OBJ=$(SRC:.f90=.o)
SUBSRC=$(addprefix src/, $(SRC))
SUBOBJ=$(addprefix obj/, $(OBJ))

ifndef dim
$(warning dimention is not set, `dim=1` will be used)
dim=1
endif
ifndef kt
$(warning kernel type is not set, `kt=n2w` will be used)
kt=n2w
endif
ifndef hf
$(warning hfac is not set, `hf=1.` will be used)
hf=1.
endif
ifndef dtype
$(warning dtype is not set, `dtype=diff` will be used)
dtype=diff
endif

%.o : src/%.f90
	@mkdir -p obj
	$(FC) $(FFLAGS) -J mod/ -o obj/$@ -c $<
IC/%.o : src/IC/%.f90
	@mkdir -p obj/IC
	@mkdir -p mod
	$(FC) $(FFLAGS) -J mod/ -o obj/$@ -c $<
kernel/%.o : src/kernel/%.f90
	@mkdir -p obj/kernel
	$(FC) $(FFLAGS) -J mod/ -o obj/$@ -c $<

execute: $(OBJ)
	@mkdir -p output
	$(FC) $(FFLAGS) -I mod/ -o $@ $(SUBOBJ)

hydroshock1:
	make && time ./execute $(dim) hydroshock 0.1 runresult.info $(kt) 1. $(dtype)

infslb:
	make && time ./execute $(dim) infslb 0.1 runresult.info $(kt) 5 $(dtype)

hc-sinx:
	make && time ./execute $(dim) hc-sinx 0.1 runresult.info $(kt) .5 $(dtype)

chi-laplace:
	make && time ./execute $(dim) chi-laplace 0.002 runresult.info $(kt) -1 $(hf) $(dtype)

pheva:
	make && time ./execute $(dim) pheva 0.002 runresult.info $(kt) -1.5 $(hf) $(dtype)

diff-laplace:
	make && time ./execute $(dim) diff-laplace .04 runresult.info $(kt) 0.1 $(hf) $(dtype)

diff-graddiv:
	make && time ./execute $(dim) diff-graddiv .04 runresult.info $(kt) 100 $(hf) $(dtype)

.PHONY: execute

clean:
	rm -rf obj/*; rm -rf mod/*; rm -r *.png; rm execute; rm -f output/*
