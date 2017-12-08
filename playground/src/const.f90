module const
  implicit none
  real, parameter :: pi = 4.*atan(1.)
  real, parameter :: eps0 = epsilon(0.)
  integer, parameter :: &
    e_none = -1,&
    ett_sin3 = 1,&
    ett_mti = 2,&
    ett_shock12 = 3,&
    ett_pulse = 4,&
    ett_ring = 5,&
    ett_soundwave = 6,&
    ett_hydroshock = 7,&
    ett_alfvenwave = 8,&
    ebc_all = 9,&
    ebc_x = 10,&
    ebc_y = 11,&
    ebc_z = 12,&
    ebc_x1 = 13,&
    ebc_y1 = 14,&
    ebc_z1 = 15,&
    ebc_x2 = 16,&
    ebc_y2 = 17,&
    ebc_z2 = 18
end module
