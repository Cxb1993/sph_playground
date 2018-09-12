module const
  implicit none
  real, parameter :: pi = 4.*atan(1.)
  real, parameter :: eps0 = epsilon(0.)

  integer, parameter :: &
    e_none = -1

  ! store of state
  integer, parameter :: &
    ec_dim = 1,&
    ec_eqs = 2,&
    ec_ddw = 3,&
    ec_ics = 4,&
    ec_adden = 5,&
    ec_artts = 6,&
    ec_coordsys = 7,&
    ec_disotropic = 8,&
    ec_dcondconst = 9,&
    ec_muzero = 10,&
    ec_tfinish = 11,&
    ec_npics = 12,&
    ec_hfac = 13,&
    ec_silent = 14,&
    ec_resolution = 15,&
    ec_spacing = 16,&
    ec_xmin = 17,&
    ec_xmax = 18,&
    ec_ymin = 19,&
    ec_ymax = 20,&
    ec_zmin = 21,&
    ec_zmax = 22,&
    ec_bordsize = 23,&
    ec_padding = 24,&
    ec_realpn = 25,&
    ec_fixedpn = 26,&
    ec_gamma = 27,&
    ec_lastprint = 28,&
    ec_usedumps = 29,&
    ec_dtprint = 30,&
    ec_process = 31,&
    ec_total = 31

  ! store of data
  integer, parameter :: &
    es_rx       = 1,&
    es_ry       = 2,&
    es_rz       = 3,&
    es_vx       = 4,&
    es_vy       = 5,&
    es_vz       = 6,&
    es_ax       = 7,&
    es_ay       = 8,&
    es_az       = 9,&
    es_m        = 10,&
    es_den      = 11,&
    es_h        = 12,&
    es_dh       = 13,&
    es_p        = 14,&
    es_u        = 15,&
    es_du       = 16,&
    es_t        = 17,&
    es_kappa    = 18,&
    es_fluxlim  = 19,&
    es_dtdx     = 20,&
    es_dtdy     = 21,&
    es_dtdz     = 22,&
    es_ddt      = 23,&
    es_bx       = 24,&
    es_by       = 25,&
    es_bz       = 26,&
    es_dbx      = 27,&
    es_dby      = 28,&
    es_dbz      = 29,&
    es_c        = 30,&
    es_om       = 31,&
    es_type     = 32,&
    es_total    = 32,&
    es_postr    = 4

  ! particles types
  integer, parameter :: &
    ept_empty     = 0,&
    ept_real      = 1,&
    ept_fixed     = 2,&
    ept_periodic  = 3

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
    ett_OTvortex    = 108,&
    ett_boilingtank = 109,&
    ett_mtilowres   = 110,&
    ett_fld_gauss   = 111

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
    eeq_magnetohydrodiffusion = 303,&
    eeq_hydrodiffusion        = 304,&
    eeq_kd2                   = 305,&
    eeq_hyrad                 = 306

  ! second derivatives
  integer, parameter :: &
    esd_fab     = 400,&
    esd_n2w     = 401,&
    esd_fw      = 402,&
    esd_2nw_ds  = 403,&
    esd_2nw_sd  = 404

  integer, parameter :: &
    edi_iso       = 500,&
    edi_aniso     = 501,&
    edi_mtih2017  = 502

  integer, parameter :: &
    epc_relaxation        = 600,&
    epc_backcompatibility = 601

  integer, parameter :: &
    eqs_hydro   = 1,&
    eqs_magneto = 2,&
    eqs_diff    = 3,&
    eqs_fld     = 4,&
    eqs_total   = 4

  character(len=*), parameter :: &
    blockFormatFlt  = "(A,A50,F16.10)",&
    blockFormatInt  = "(A,A50,I8)",&
    blockFormatStr  = "(A,A50,A)",&
    blockFormatStr2 = "(A,A50,A)"

  integer, parameter :: &
    enl_l1 = 1,&
    enl_l2 = 2

end module
