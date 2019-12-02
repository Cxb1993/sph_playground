# SPH playground
Fortran code and python runners for SPH

This code was written by me during PhD at MoCA.

The SPH core
---
[sph_playground/playground](https://github.com/Evedel/sph_playground/tree/master/playground)

The code meant to be self-explanatory.

Clearly, it is not.

All-mighty `Fortran` SPH code.
---
 [playground/src](https://github.com/Evedel/sph_playground/tree/master/playground/src)

`Python` drivers for SPH engine.
---
[playground/runner](https://github.com/Evedel/sph_playground/tree/master/playground/runner)

The idea was
---
To write readable python-like-almost-human-like code in runners.

All runners are based on [playground/runner/fortdriver](https://github.com/Evedel/sph_playground/tree/master/playground/runner/fortdriver).

1. Fortdriver reads all the constants and enumerators from [const.f90](https://github.com/Evedel/sph_playground/blob/master/playground/src/const.f90)
2. Fortdriver compiles SPH engine
3. Fortdriver run the initial setup of SPH engine with respect to `Context.setup`
4. Fortdriver may access all the results of the run and modify dumps as you like
5. Additional runs after dumps were modified? Sure thing.

Typical examples are [convergence of aisotropic diffusion on a glass like latice for different methods and initial conditions](https://github.com/Evedel/sph_playground/blob/master/playground/runner/000_anisopaper_v1_gls_iso.py) ([paper.pdf](https://arxiv.org/pdf/1812.04006.pdf))
and [flux-limited diffusion radiation on a shock tube](https://github.com/Evedel/sph_playground/blob/master/playground/runner/004.3_fluxlimdif_shock.py).

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

What is in this code
---
- N-deminsional SPH in two loops
- bunch of kernels that can be chosen on compile time only (FortDriver can handle this) ([playground/src/kernel](https://github.com/Evedel/sph_playground/tree/master/playground/src/kernel))
- bunch of second derivative methods that can be chosen on run time and implemented with the function pointers ([playground/src/kernel/kernel.f90](https://github.com/Evedel/sph_playground/blob/master/playground/src/kernel/kernel.f90))
- close-packed, glass-like and, of course, uniform initial latices ([playground/src/IC](https://github.com/Evedel/sph_playground/tree/master/playground/src/IC))
- bunch of specific initial conditions, which are deprecated as it is the age of FortDriver now ([playground/src/IC/IC.f90](https://github.com/Evedel/sph_playground/blob/master/playground/src/IC/IC.f90))
- own implementations of [appendable arrays](https://github.com/Evedel/sph_playground/blob/master/playground/src/utils/arrayresize.f90), 
[(python-like) lists](https://github.com/Evedel/sph_playground/blob/master/playground/src/utils/list.f90), 
[dictionaries](https://github.com/Evedel/sph_playground/blob/master/playground/src/utils/map.f90) and [others](https://github.com/Evedel/sph_playground/tree/master/playground/src/utils)
- free, periodic, mirroring and fixed border conditions ([playground/src/BC.f90](https://github.com/Evedel/sph_playground/blob/master/playground/src/BC.f90))
- hydrodynamics, magnetohydrodynamics, diffusion, radiative hydrodynamics ([playground/src/circuit2.f90](https://github.com/Evedel/sph_playground/blob/master/playground/src/circuit2.f90))
- statefull hdf5 dumps and the restore of a stopped/modified simmulation ([playground/src/state.f90](https://github.com/Evedel/sph_playground/blob/master/playground/src/state.f90))
- Leapfrog and Super-Time-Stepping integrators ([playground/src/sts_integrator.f90](https://github.com/Evedel/sph_playground/blob/master/playground/src/sts_integrator.f90))
- KDTree based neighbour search ([playground/src/utils/neighboursearch.f90](https://github.com/Evedel/sph_playground/blob/master/playground/src/utils/neighboursearch.f90))
- cartesian and cilindrical coordinate systems

What to do if you want to use this code
---
JUST DON'T
---

The code was implemented only for research purposes.

There are plenty more fish in the sea.

Shell scripts for kernels and derivatives tests
---
[sph_playground/scripts](https://github.com/Evedel/sph_playground/tree/master/scripts)

M_{6/2} kernel for [Phantom](https://bitbucket.org/danielprice/phantom/src/master/)
---
[sph_playground/kernel_zipm6.f90](https://github.com/Evedel/sph_playground/blob/master/kernel_zipm6.f90)

Division vs multiplication test
---
A simple test which checks that in `Fortran` it is indeed up to `30%` faster to use `*0.25` rather than `\4.`

[sph_playground/test_4vs025.f90](https://github.com/Evedel/sph_playground/blob/master/test_4vs025.f90)

Compiled with  `gfortran -O4 -Wall test_4vs025.f90 -fopenmp -o test_4vs025`

Memory allignment test
---
[sph_playground/test_soa_aos.f90](https://github.com/Evedel/sph_playground/blob/master/test_soa_aos.f90)

Compiled with  `gfortran -O4 -Wall test_soa_aos.f90 -fopenmp -o test_soa_aos`

<img src="https://github.com/Evedel/sph_playground/blob/master/mem.png" width="350" height="400">

This test is simple, which is why the code is just a horrible boilerplate.

The idea of the test is to determine which memory allignment is the fastest for SPH in Fortran.

We iterate over the array of "particles" with the number of particles N=10^4...10^12.

For each particle, we store mass, force, acceleration, velocity, and position in 3D. All are randomly generated numbers outside of the time measure regions.

For each particle we do some "physics" *(by the way this is the example of the winner -- AOS)*

```
P(i).Acceleration.x = P(i).Force.x / P(i).Mass
P(i).Acceleration.y = P(i).Force.y / P(i).Mass
P(i).Acceleration.z = P(i).Force.z / P(i).Mass

P(i).Velocity.x = P(i).Velocity.x + P(i).Acceleration.x * dt
P(i).Velocity.y = P(i).Velocity.y + P(i).Acceleration.y * dt
P(i).Velocity.z = P(i).Velocity.z + P(i).Acceleration.z * dt

P(i).Position.x = P(i).Position.x + P(i).Velocity.x * dt
P(i).Position.y = P(i).Position.y + P(i).Velocity.x * dt
P(i).Position.z = P(i).Position.z + P(i).Velocity.z * dt
```
The main reason for such "physics" is to enforce the interconnection of preperties and simultaneous read/write operations for each particle.

We check sequential `i=1...N` and random accesses. For random access, we create a random array of indexes `random(indexes[1...N])` that need to be updated during the main loop. After that, we iterate over the indexes in the array. The number of operations is `O(N)` again.

The memory alignments are:

SOA (1 loop) -- Structure Of Arrays.
```
Storage {
  mass(1...N)
  force_x(1...N)
  ...
  position_z(1...N)
}
```
SOA (3 loops) -- same as SOA (1 loop), but we used 3 smaller main loops in particle update. Because why not.

AOS -- Array Of Structures.
```
Storage(1...N) = [
  Particle {
    mass
    force {
      x
      y
      z
    }
    ...
    position {
      ...
    }
  }
]
```
NDA -- N-N-Dimensional Arrays
```
Acceleration(x...z, 1...N) = Force(x...z, 1...N)/Mass(1...N)
Velocity(x...z, 1...N) = Velocity(x...z, 1...N) + Acceleration(x...z, 1...N)*dt
Position(x...z, 1...N) = Position(x...z, 1...N) + Velocity(x...z, 1...N)*dt
```
SAR -- Super-Array, where the second index is the number of a particle and the first index is the number of a property (mass, force_x, ... position_z)
```
Storage(1...13, 1...N)
```

#### Even though AOS is believed to be the slowest memory alignment, it is the fastest in practice for `Fortran` codes (especially for random access, which is dominant in SPH). Additionally, AOS is the same is SAR in `Fortran`, without computational overhead but much nicer looking (see the example of the "physics").
