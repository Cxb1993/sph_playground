The code meant to be self-explanatory.

Clearly, it is not.

All-mighty `Fortran` SPH code.
---
 [playground/src](https://github.com/Evedel/sph_playgroung/tree/master/playground/src)

`Python` drivers for SPH engine.
---
[playground/runner](https://github.com/Evedel/sph_playgroung/tree/master/playground/runner)

The idea was
---
To write readable python-like-almost-human-like code in runners.

All runners are based on [playground/runner/fortdriver](https://github.com/Evedel/sph_playgroung/tree/master/playground/runner/fortdriver).

1. Fortdriver reads all the constants and enumerators from [const.f90](https://github.com/Evedel/sph_playgroung/blob/master/playground/src/const.f90)
2. Fortdriver compiles SPH engine
3. Fortdriver run the initial setup of SPH engine with respect to `Context.setup`
4. Fortdriver may access all the results of the run and modify dumps as you like
5. Additional runs after dumps were modified? Sure thing.

Typical examples are [convergence of aisotropic diffusion on a glass like latice for different methods and initial conditions](https://github.com/Evedel/sph_playgroung/blob/master/playground/runner/000_anisopaper_v1_gls_iso.py) ([paper.pdf](https://arxiv.org/pdf/1812.04006.pdf))
and [flux-limited diffusion radiation on a shock tube](https://github.com/Evedel/sph_playgroung/blob/master/playground/runner/004.3_fluxlimdif_shock.py).

It is possible to run the code without FortDriver in backward compatibility mode, but it requires a lot of hardcoded lines and ambiguous execute commands, for example
```
>> make kernel=m6 useomp=t debug=f

>> ./execute --xmax 1.0 --nsnapshots 10.0 --silent no --au 1.0
            --equations hydro --influencefile influence.info
            --xmin -1.0 --artts yes --adden yes --resolution 512
            --hfac 1.0 --coordsys cartesian --ics shock12
            --resultfile result.info --process relaxation
            --usedumps yes --dim 2.0 --tfinish 20.0 --ddw fab
```
What is in code
- N-deminsional SPH in two loops
- bunch of kernels that can be chosen on compile time only (FortDriver can handle this) ([playground/src/kernel](https://github.com/Evedel/sph_playgroung/tree/master/playground/src/kernel))
- bunch of second derivative methods that can be chosen on run time and implemented with the function pointers ([playground/src/kernel/kernel.f90](https://github.com/Evedel/sph_playgroung/blob/master/playground/src/kernel/kernel.f90))
- close-packed, glass-like and, of course, uniform initial latices ([playground/src/IC](https://github.com/Evedel/sph_playgroung/tree/master/playground/src/IC))
- bunch of specific initial conditions, which are deprecated as it is the age of FortDriver now ([playground/src/IC/IC.f90](https://github.com/Evedel/sph_playgroung/blob/master/playground/src/IC/IC.f90))
- own implementations of [appendable arrays](https://github.com/Evedel/sph_playgroung/blob/master/playground/src/utils/arrayresize.f90), 
[(python-like) lists](https://github.com/Evedel/sph_playgroung/blob/master/playground/src/utils/list.f90), 
[dictionaries](https://github.com/Evedel/sph_playgroung/blob/master/playground/src/utils/map.f90) and [others](https://github.com/Evedel/sph_playgroung/tree/master/playground/src/utils)
- free, periodic, mirroring and fixed border conditions ([playground/src/BC.f90](https://github.com/Evedel/sph_playgroung/blob/master/playground/src/BC.f90))
- hydrodynamics, magnetohydrodynamics, diffusion, radiative hydrodynamics ([playground/src/circuit2.f90](https://github.com/Evedel/sph_playgroung/blob/master/playground/src/circuit2.f90))
- statefull hdf5 dumps and the restore of a stopped/modified simmulation ([playground/src/state.f90](https://github.com/Evedel/sph_playgroung/blob/master/playground/src/state.f90))
- Leapfrog and Super-Time-Stepping integrators ([playground/src/sts_integrator.f90](https://github.com/Evedel/sph_playgroung/blob/master/playground/src/sts_integrator.f90))
- KDTree based neighbour search ([playground/src/utils/neighboursearch.f90](https://github.com/Evedel/sph_playgroung/blob/master/playground/src/utils/neighboursearch.f90))
- cartesian and cilindrical coordinate systems

What to do if you want to use this code
---
JUST DON'T
---

The code was implemented only for research purposes.

There are plenty more fish in the sea.
