module const
  implicit none
  real, parameter :: pi = 4.*atan(1.)
  real, parameter :: eps0 = epsilon(0.)

  integer, parameter :: &
    e_none = -1

  ! initial setups
  integer, parameter :: &
    ett_sin3 = 1,&
    ett_mti = 2,&
    ett_shock12 = 3,&
    ett_pulse = 4,&
    ett_ring = 5,&
    ett_soundwave = 6,&
    ett_hydroshock = 7,&
    ett_alfvenwave = 8,&
    ett_OTvortex = 9

  ! borders definitions
  integer, parameter :: &
    ebc_all = 100,&
    ebc_x = 101,&
    ebc_y = 102,&
    ebc_z = 103,&
    ebc_x1 = 104,&
    ebc_y1 = 105,&
    ebc_z1 = 106,&
    ebc_x2 = 107,&
    ebc_y2 = 108,&
    ebc_z2 = 109

  ! set of equations 
  integer, parameter :: &
    eeq_hydro = 201,&
    eeq_diffusion = 202,&
    eeq_magnetohydro = 203,&
    eeq_magnetohydrodiffusion = 204
end module
