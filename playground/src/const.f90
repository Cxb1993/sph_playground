module const
  implicit none
  real, parameter :: pi = 4.*atan(1.)
  real, parameter :: eps0 = epsilon(0.)

  integer, parameter :: &
    e_none = -1

  ! store of data
  integer, parameter :: &
    es_posx   = 1,&
    es_posy   = 2,&
    es_posz   = 3,&
    es_velx   = 4,&
    es_vely   = 5,&
    es_velz   = 6,&
    es_accx   = 7,&
    es_accy   = 8,&
    es_accz   = 9,&
    es_total  = 9

  ! initial setups
  integer, parameter :: &
    ett_sin3        = 100,&
    ett_mti         = 101,&
    ett_shock12     = 102,&
    ett_pulse       = 103,&
    ett_ring        = 104,&
    ett_soundwave   = 105,&
    ett_hydroshock  = 106,&
    ett_alfvenwave  = 107,&
    ett_OTvortex    = 108

  ! borders definitions
  integer, parameter :: &
    ebc_all = 200,&
    ebc_x   = 201,&
    ebc_y   = 202,&
    ebc_z   = 203,&
    ebc_x1  = 204,&
    ebc_y1  = 205,&
    ebc_z1  = 206,&
    ebc_x2  = 207,&
    ebc_y2  = 208,&
    ebc_z2  = 209

  ! set of equations
  integer, parameter :: &
    eeq_hydro                 = 300,&
    eeq_diffusion             = 301,&
    eeq_magnetohydro          = 302,&
    eeq_magnetohydrodiffusion = 303

  integer, parameter :: &
    esd_fab = 400,&
    esd_n2w = 401,&
    esd_fw  = 402,&
    esd_2nw = 403

  character(len=*), parameter :: &
    blockFormatFlt  = "(A,A35,F10.7)",&
    blockFormatInt  = "(A,A34,I6)",&
    blockFormatStr  = "(A,A35,A)",&
    blockFormatStr2 = "(A,A31,A)"

end module
