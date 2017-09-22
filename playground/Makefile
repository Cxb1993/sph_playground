#! /usr/bin/make -f

ifndef kernel
$(warning base kernel name not provided, `zipm6` will be used)
kernel=zipm6
endif

SRC= \
	utils/arrayresize.f90\
	utils/timing.f90\
	utils/map.f90\
	const.f90\
	state.f90

SRC+=kernel/base_$(kernel).f90

SRC+=kernel/kernel.f90\
	utils/neighboursearch.f90\
	BC.f90\
	printer.f90\
	err_shocktube.f90\
	errteylor.f90\
	args.f90\
	eos.f90\
	IC/initpositions.f90\
	IC/IC.f90\
	IC/initvariables.f90\
	circuit1.f90\
	circuit2.f90\
	errcalc.f90\
	iterator.f90\
	main.f90

OBJ=$(SRC:.f90=.o)
SUBSRC=$(addprefix src/, $(SRC))
SUBOBJ=$(addprefix obj/, $(OBJ))

# $(warning $(SRC))

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
ifndef taskcheck
taskcheck=chi-laplace
endif
ifndef debug
debug=no
endif
ifndef useomp
useomp=f
endif

FC=gfortran

FFLAGS=-fdefault-real-8 -fdefault-double-8\
 				-O4

ifeq ($(debug),yes)
FFLAGS=-fdefault-real-8 -fdefault-double-8\
 				-O0 -g \
				-Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace\
				-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow,invalid -finit-real=nan
endif

# ifeq ($(useomp),t)
# 	FFLAGS+=-fopenmp
# endif

FFLAGS+=-fopenmp

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

utils/%.o : src/utils/%.f90
	@mkdir -p obj/utils
	$(FC) $(FFLAGS) $(modflagout) mod/ -o obj/$@ -c $<

execute: $(OBJ)
	@mkdir -p output
	$(FC) $(FFLAGS) $(modflaginp) mod/ -o $@ $(SUBOBJ)

memleakcheck:
	# --suppressions=.configure \
	make && valgrind -v --track-origins=yes\
			--leak-check=full --show-leak-kinds=all \
			--gen-suppressions=all \
			./execute --dim $(dim)\
			--tasktype $(taskcheck) --spacing .01\
			--errfilename runresult.info --kerneltype $(kt)\
			--tfinish .2 --hfac $(hf) --difftype $(dtype)\
			--silent $(silent)

hydroshock:
	make debug=$(debug) && time ./execute --dim $(dim)\
			--tasktype hydroshock --spacing 0.0005\
			--errfilename runresult.info --kerneltype $(kt)\
			--tfinish .2 --hfac $(hf) --difftype $(dtype)\
			--silent $(silent)

soundwave:
	make debug=$(debug) && time ./execute --dim $(dim)\
			--tasktype soundwave --spacing .1\
			--errfilename runresult.info --kerneltype $(kt)\
			--tfinish .2 --hfac $(hf) --difftype $(dtype)\
			--silent $(silent)

infslb:
	make debug=$(debug) && time ./execute $(dim) infslb 0.1 runresult.info $(kt) 5 $(dtype) $(silent)

hc-sinx:
	make debug=$(debug) && time ./execute --dim $(dim)\
			--tasktype hc-sinx --spacing .1\
			--errfilename runresult.info --kerneltype $(kt)\
			--tfinish 5. --hfac $(hf) --difftype $(dtype)\
			--silent $(silent)

pheva:
	make debug=$(debug) && time ./execute $(dim) pheva 0.002 runresult.info $(kt) -1.5 $(hf) $(dtype) $(silent)

diff-laplace:
	make debug=$(debug) && time ./execute --dim $(dim)\
			--tasktype diff-laplace --spacing .06\
			--errfilename runresult.info --kerneltype $(kt)\
			--tfinish 100 --hfac $(hf) --difftype $(dtype)\
			--silent $(silent)

diff-graddiv:
	make debug=$(debug) && time ./execute $(dim) diff-graddiv .05 runresult.info $(kt) 100 $(hf) $(dtype) $(silent)

chi-laplace:
	make debug=$(debug) && time ./execute --dim $(dim)\
			--tasktype chi-laplace --spacing .06\
			--errfilename runresult.info --kerneltype $(kt)\
			--tfinish -1 --hfac $(hf) --difftype $(dtype)\
			--silent $(silent)

chi-graddiv:
	make debug=$(debug) && time ./execute --dim $(dim)\
			--tasktype chi-graddiv --spacing .06\
			--errfilename runresult.info --kerneltype $(kt)\
			--tfinish -1 --hfac $(hf) --difftype $(dtype)\
			--silent $(silent)

inflience-graddiv:
	make debug=$(debug) && time ./execute --dim $(dim) \
			--tasktype chi-graddiv --spacing 0.06 \
			--errfilename runresult.info --kerneltype $(kt)\
			--tfinish 100 --hfac $(hf) --difftype $(dtype) \
			--kerninfluencefile influence.info --silent $(silent)

inflience-laplace:
	make debug=$(debug) && time ./execute --dim $(dim) \
			--tasktype chi-laplace --spacing 0.06 \
			--errfilename runresult.info --kerneltype $(kt)\
			--tfinish 100 --hfac $(hf) --difftype $(dtype) \
			--kerninfluencefile output/influence.info --silent $(silent)

.PHONY: execute

clean:
	rm -rf obj/*; rm -rf mod/*; rm -r *.png; rm execute; rm -f output/*
