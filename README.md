# SPH playgroung
fortran code and python runners for SPH

This code was written and used during my PhD at MoCA.

---
The SPH core

[sph_playgroung/playground](https://github.com/Evedel/sph_playgroung/tree/master/playground)

---
Shell scripts for kernels and derivatives tests

[sph_playgroung/scripts](https://github.com/Evedel/sph_playgroung/tree/master/scripts)

---
M_{6/2} kernel for [Phantom](https://bitbucket.org/danielprice/phantom/src/master/)

[sph_playgroung/kernel_zipm6.f90](https://github.com/Evedel/sph_playgroung/blob/master/kernel_zipm6.f90)

---
Division vs multiplicaion test

A simple test which checks that in `fortran` it is indeed upto `30%` faster to use `*0.25` rather than `\4.`

[sph_playgroung/test_4vs025.f90](https://github.com/Evedel/sph_playgroung/blob/master/test_4vs025.f90)

---
Memory allignment test

[sph_playgroung/test_soa_aos.f90](https://github.com/Evedel/sph_playgroung/blob/master/test_soa_aos.f90)

---
In its current form the code supposed to be run with python scripts
