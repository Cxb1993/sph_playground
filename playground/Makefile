#! /usr/bin/make -f

ifndef kernel
$(warning base kernel name not provided, `zipm6` will be used)
kernel=zipm6
endif

SRC= \
	utils/arrayresize.f90\
	utils/timing.f90\
	utils/map.f90\
	utils/list.f90\
	utils/kdtree2.f90\
	const.f90\
	utils/tables.f90\
	utils/stretchmap.f90\
	state.f90\
	BC.f90\
	preruncheck.f90

SRC+=kernel/base_$(kernel).f90

SRC+=kernel/kernel.f90\
	utils/neighboursearch.f90\
	printer.f90\
	err_shocktube.f90\
	errteylor.f90\
	args.f90\
	eos.f90\
	IC/uniform.f90\
	IC/rhofuncs.f90\
	IC/setup.f90\
	IC/initvariables.f90\
	circuit1.f90\
	circuit2.f90\
	errcalc.f90\
	iterator.f90\
	utils/clfortran.f90\
	cltest.f90\
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
debug=f
endif
ifndef useomp
useomp=t
endif

ifndef xperf
xperf=f
endif

FC=gfortran

FFLAGS=-fdefault-real-8 -fdefault-double-8 \
 				-flto -Ofast

ifeq ($(debug),t)
FFLAGS=-fdefault-real-8 -fdefault-double-8\
 				-O0 -g \
				-Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace\
				-ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow,invalid -finit-real=nan
endif
ifeq ($(debug),p)
FFLAGS=-pg -no-pie -fdefault-real-8 -fdefault-double-8 -flto -Ofast
endif

ifeq ($(useomp),t)
	FFLAGS+=-fopenmp
endif

ifeq ($(xperf),f)
	FFLAGS+=-O4
endif
ifeq ($(xperf),t)
	FFLAGS+=-fprofile-generate=profiler
endif
ifeq ($(xperf),c)
	FFLAGS+=-fprofile-use=profiler -fprofile-correction
endif
# FFLAGS+=-fopenmp

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
			--errfilename result.info --kerneltype $(kt) \
			--tfinish .05 --hfac 1.2 --difftype diff --silent no

.PHONY: execute

clean:
	rm -rf obj/*; rm -rf mod/*; rm -r *.png; rm execute; rm -f output/*
