#! /usr/bin/make -f

FC=gfortran
FFLAGS=-fdefault-real-8 -fdefault-double-8\
 				-O4 \
				-Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace\
				-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow,invalid -finit-real=nan\
				-fopenmp

# FC := pgfortran
# FFLAGS=-r8 -Mallocatable=03 -mp -Mconcur=allcores -Mipa -Mprof -Minfo
				# -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow,invalid -finit-real=nan
				# -fopenmp
				# -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace
				# -g

# kernel/base_gaus.f90\
# kernel/n2movedgaus.f90\

SRC=\
	utils.f90\
	timing.f90\
	const.f90\
	kernel/base_cubic.f90\
	kernel/base_quintic.f90\
	kernel/base_external.f90\
	kernel/base_sinq.f90\
	kernel/n2ext.f90\
	kernel/kernel.f90\
	neighboursearch.f90\
	BC.f90\
	errteylor.f90\
	args.f90\
	printer.f90\
	eos.f90\
	IC/uniform.f90\
	IC/semiuniform.f90\
	IC/IC.f90\
	circuit1.f90\
	circuit2.f90\
	errcalc.f90\
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
ifndef silent
$(warning silent is not set, `silent=no` will be used)
silent=yes
endif

ifeq ($(FC),gfortran)
modflagout=-J
modflaginp=-I
endif
ifeq ($(FC),pgfortran)
modflagout=-module
modflaginp=-module
endif

%.o : src/%.f90
	@mkdir -p obj
	$(FC) $(FFLAGS) $(modflagout) mod/ -o obj/$@ -c $<

IC/%.o : src/IC/%.f90
	@mkdir -p obj/IC
	@mkdir -p mod
	$(FC) $(FFLAGS) $(modflagout) mod/ -o obj/$@ -c $<

kernel/%.o : src/kernel/%.f90
	@mkdir -p obj/kernel
	$(FC) $(FFLAGS) $(modflagout) mod/ -o obj/$@ -c $<

execute: $(OBJ)
	@mkdir -p output
	$(FC) $(FFLAGS) $(modflaginp) mod/ -o $@ $(SUBOBJ)

hydroshock1:
	make && time ./execute $(dim) hydroshock 0.1 runresult.info $(kt) 1. $(dtype) $(silent)

infslb:
	make && time ./execute $(dim) infslb 0.1 runresult.info $(kt) 5 $(dtype) $(silent)

hc-sinx:
	make && time ./execute $(dim) hc-sinx 0.1 runresult.info $(kt) .5 $(dtype) $(silent)

pheva:
	make && time ./execute $(dim) pheva 0.002 runresult.info $(kt) -1.5 $(hf) $(dtype) $(silent)

diff-laplace:
	make && time ./execute $(dim) diff-laplace .04 runresult.info $(kt) 0.1 $(hf) $(dtype) $(silent)

diff-graddiv:
	make && time ./execute $(dim) diff-graddiv .05 runresult.info $(kt) 100 $(hf) $(dtype) $(silent)

chi-laplace:
	make && time ./execute $(dim) chi-laplace 0.06 runresult.info $(kt) -1 $(hf) $(dtype) $(silent)

chi-graddiv:
	make && time ./execute $(dim) chi-graddiv 0.06 runresult.info $(kt) -1 $(hf) $(dtype) $(silent)

inflience:
	make && time ./execute --dim $(dim) \
			--tasktype chi-laplace --spacing 0.06 \
			--errfilename runresult.info --kerneltype $(kt)\
			--tfinish 100 --hfac $(hf) --difftype $(dtype) \
			--kerninfluencefile f.info --silent $(silent)

.PHONY: execute

clean:
	rm -rf obj/*; rm -rf mod/*; rm -r *.png; rm execute; rm -f output/*
