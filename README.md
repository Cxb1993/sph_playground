# SPH playground
Fortran code and python runners for SPH

This code was written by me during PhD at MoCA.

The SPH core
---
For more info go to [sph_playgroung/playground](https://github.com/Evedel/sph_playgroung/tree/master/playground)

Shell scripts for kernels and derivatives tests
---
[sph_playgroung/scripts](https://github.com/Evedel/sph_playgroung/tree/master/scripts)

M_{6/2} kernel for [Phantom](https://bitbucket.org/danielprice/phantom/src/master/)
---
[sph_playgroung/kernel_zipm6.f90](https://github.com/Evedel/sph_playgroung/blob/master/kernel_zipm6.f90)

Division vs multiplication test
---
A simple test which checks that in `Fortran` it is indeed up to `30%` faster to use `*0.25` rather than `\4.`

[sph_playgroung/test_4vs025.f90](https://github.com/Evedel/sph_playgroung/blob/master/test_4vs025.f90)
Compiled with  `gfortran -O4 -Wall test_4vs025.f90 -fopenmp -o test_4vs025`

Memory allignment test
---
[sph_playgroung/test_soa_aos.f90](https://github.com/Evedel/sph_playgroung/blob/master/test_soa_aos.f90)

Compiled with  `gfortran -O4 -Wall test_soa_aos.f90 -fopenmp -o test_soa_aos`

<img src="https://github.com/Evedel/sph_playgroung/blob/master/mem.png" width="350" height="400">

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
