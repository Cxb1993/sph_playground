module state
  use omp_lib

  use const
  use errprinter, only: error, warning

  implicit none

  public :: set_equations, get_equations, &
            setddwtype, getddwtype, getdim, setdim, &
            ginitvar, sinitvar, setAdvancedDensity, getAdvancedDensity,&
            setArtificialTerms, getArtificialTerms, &
            scoordsys, gcoordsys, &
            setdiffconductivity, getdiffconductivity, &
            setdiffisotropic, getdiffisotropic, &
            setmhdmagneticpressure, getmhdmagneticpressure, &
            printstate, setnpics, settfinish, gethfac, gettfinish,&
            getnpics, sethfac, setsilentmode, getsilentmode,&
            getspacing, setspacing, setresolution, getresolution,&
            setBorders, getBorders,&
            setPartNumber, getPartNumber,&
            setState, getState,&
            setGamma, getGamma,&
            setLastPrint,getLastPrint,&
            setUseDumps, getUseDumps,&
            setdtprint, getdtprint,&
            getEqComponent,&
            setProcess, getProcess,&
            setArtTermCond, getArtTermCond,&
            clearState,&
            setStateVal, getStateVal,&
            splacement

  private
  save
    real :: statevars(ec_total)

    interface getStateVal
       module procedure getStateVal_r, getStateVal_i
    end interface

  contains
    subroutine clearState()
      statevars(:) = -1
    end subroutine clearState

    subroutine setStateVal(label, val)
    integer, intent(in) :: label
    real, intent(in) :: val

      statevars(label) = val
    end subroutine setStateVal

    subroutine getStateVal_r(label, val)
    integer, intent(in) :: label
    real, intent(out) :: val
      val = statevars(label)
    end subroutine getStateVal_r

    subroutine getStateVal_i(label, val)
    integer, intent(in) :: label
    integer, intent(out) :: val
      val = int(statevars(label))
      if ((label /= ec_dim).and.&
          (label /= ec_lastprint).and.&
          (label /= ec_fixedpn).and.&
          (label /= ec_realpn).and.&
          (label /= ec_usedumps).and.&
          (label /= ec_ststype).and.&
          (label /= ec_stsfixeds).and.&
          (label /= ec_silent)) then
          call warning("Integer conversion were used unexpectedly", label, __FILE__, __LINE__)
          call warning("Value before", statevars(label), __FILE__, __LINE__)
          call warning("Value after", int(statevars(label)), __FILE__, __LINE__)
      end if
    end subroutine getStateVal_i

    subroutine setdim(d)
      integer, intent(in) :: d
      if ((d < 1).or.(d > 3)) then
        call error('Dimmention is wrong', d, __FILE__, __LINE__)
      else
        statevars(ec_dim) = d
      end if
    end subroutine
    pure subroutine getdim(d)
      integer, intent(out) :: d
      d = int(statevars(ec_dim))
    end subroutine

    subroutine set_equations(itt)
      character (len=*), intent(in) :: itt
      select case(itt)
      case('hydro')
        statevars(ec_eqs) = eeq_hydro
      case('magnetohydro')
        statevars(ec_eqs) = eeq_magnetohydro
      case('diffusion')
        statevars(ec_eqs) = eeq_diffusion
      case('mhd')
        statevars(ec_eqs) = eeq_magnetohydrodiffusion
      case('hydi')
        statevars(ec_eqs) = eeq_hydrodiffusion
      case('kd2')
        statevars(ec_eqs) = eeq_kd2
      case('hyrad')
        statevars(ec_eqs) = eeq_hyrad
      case default
        statevars(ec_eqs) = e_none
        call error('Equation is not set', itt, __FILE__, __LINE__)
      end select
    end subroutine set_equations
    pure subroutine get_equations(ott)
      integer, intent(out) :: ott
      ott = int(statevars(ec_eqs))
    end subroutine get_equations
    subroutine getEqComponent(eqSet)
      integer, intent(out) :: eqSet(eqs_total)
      eqSet(:) = 0

      select case(int(statevars(ec_eqs)))
      case(eeq_hydro)
        eqSet(eqs_hydro) = 1
      case(eeq_magnetohydro)
        eqSet(eqs_hydro)   = 1
        eqSet(eqs_magneto) = 1
      case(eeq_diffusion)
        eqSet(eqs_diff) = 1
      case(eeq_hyrad)
        eqSet(eqs_fld)     = 1
        eqSet(eqs_hydro)   = 1
        eqSet(eqs_radexch) = 1
      case(eeq_magnetohydrodiffusion)
        eqSet(eqs_hydro)   = 1
        eqSet(eqs_magneto) = 1
        eqSet(eqs_diff)    = 1
      case(eeq_hydrodiffusion)
        eqSet(eqs_hydro) = 1
        eqSet(eqs_diff)  = 1
      case(eeq_manual)
        if (int(statevars(ec_eqondiff))==1)     eqSet(eqs_diff)    = 1
        if (int(statevars(ec_eqonhydro))==1)    eqSet(eqs_hydro)   = 1
        if (int(statevars(ec_eqonfluxlim))==1)  eqSet(eqs_fld)     = 1
        if (int(statevars(ec_eqonradexch))==1)  eqSet(eqs_radexch) = 1
        if (int(statevars(ec_eqonsts))==1)      eqSet(eqs_sts)     = 1
      case default
        call warning("Equations are not set. Returning all zeros.", int(statevars(ec_eqs)), __FILE__, __LINE__)
      end select
    end subroutine getEqComponent

    subroutine setddwtype(itt)
      character (len=*), intent(in) :: itt
      select case(itt)
      case('n2w')
        statevars(ec_ddw) = esd_n2w
      case('fab')
        statevars(ec_ddw) = esd_fab
      case('2nw-+')
        statevars(ec_ddw) = esd_2nw_ds
      case('2nw+-')
        statevars(ec_ddw) = esd_2nw_sd
      case('fw')
        statevars(ec_ddw) = esd_fw
      case default
        statevars(ec_ddw) = e_none
        call error("Second derivative kernel is not set", itt, __FILE__, __LINE__)
      end select
    end subroutine setddwtype
    pure subroutine getddwtype(ott)
      integer, intent(out) :: ott
      ott = int(statevars(ec_ddw))
    end subroutine getddwtype

    subroutine sinitvar(itt, ott)
      character (len=*), intent(in) :: itt
      real, intent(out) :: ott

      select case(itt)
      case('sin3')
        ott = ett_sin3
      case('mti')
        ott = ett_mti
      case('mtilowres')
        ott = ett_mtilowres
      case('shock12')
        ott = ett_shock12
      case('pulse')
        ott = ett_pulse
      case('ring')
        ott = ett_ring
      case('soundwave')
        ott = ett_soundwave
      case('hydroshock')
        ott = ett_hydroshock
      case('alfvenwave')
        ott = ett_alfvenwave
      case('otvortex')
        ott = ett_OTvortex
      case('boilingtank')
        ott = ett_boilingtank
      case('fld_gauss')
        ott = ett_fld_gauss
      case('none')
        ott = e_none
      case default
        ott = e_none
        call error('There is no such initial variable', itt, __FILE__, __LINE__)
      end select

      if (ott /= e_none) then
        call warning("The ICS and setupV2 is depricated.", "Use place + proces methods.", __FILE__, __LINE__)
      end if
    end subroutine
    pure subroutine ginitvar(ott)
     integer, intent(out) :: ott
     ott = int(statevars(ec_ics))
   end subroutine

   subroutine setAdvancedDensity(iad)
     character (len=*), intent(in) :: iad
     select case(iad)
     case('yes')
       statevars(ec_adden) = eif_yes
     case('no')
       statevars(ec_adden) = eif_no
     case('')
       statevars(ec_adden) = eif_yes
     case default
       print *, 'There is no such advanced density choise: ', iad
       stop
     end select
   end subroutine
   pure subroutine getAdvancedDensity(oad)
     integer, intent(out) :: oad
     oad = int(statevars(ec_adden))
   end subroutine

   subroutine setArtificialTerms(iat)
     character (len=*), intent(in) :: iat
     select case(iat)
     case('yes')
       statevars(ec_artts) = eif_yes
     case('no')
       statevars(ec_artts) = eif_no
     case('')
       statevars(ec_artts) = eif_yes
     case default
       call error("There is no such default Artificail Term", iat, __FILE__, __LINE__)
     end select
   end subroutine
   pure subroutine getArtificialTerms(oat)
     integer, intent(out) :: oat
     oat = int(statevars(ec_artts))
   end subroutine

  subroutine setArtTermCond(iat)
    real, intent(in) :: iat
    statevars(ec_au) = iat
  end subroutine setArtTermCond
  pure subroutine getArtTermCond(oat)
    real, intent(out) :: oat
    oat = statevars(ec_au)
  end subroutine getArtTermCond

   subroutine scoordsys(ics)
     character (len=*), intent(in) :: ics
     select case(ics)
     case('cartesian')
       statevars(ec_coordsys) = 1
     case('cylindric')
       statevars(ec_coordsys) = 2
     case('')
       statevars(ec_coordsys) = 1
     case default
       print *, 'There is no such coordinate system: ', ics
       stop
      end select
   end subroutine
   pure subroutine gcoordsys(ocs)
     integer, intent(out) :: ocs
     ocs = int(statevars(ec_coordsys))
   end subroutine

   subroutine setdiffconductivity(i)
     real, intent(in) :: i
     statevars(ec_dcondconst) = i
   end subroutine
   pure subroutine getdiffconductivity(o)
     real, intent(out) :: o
     o = statevars(ec_dcondconst)
   end subroutine

   subroutine setdiffisotropic(i)
     character (len=*), intent(in) :: i
     select case(i)
     case('yes')
       statevars(ec_disotropic) = eif_yes
     case('no')
       statevars(ec_disotropic) = eif_no
     case('')
       statevars(ec_disotropic) = eif_yes
     case default
       call error('There is no such diso mode: ', i, __FILE__, __LINE__)
     end select
   end subroutine
   pure subroutine getdiffisotropic(o)
     integer, intent(out) :: o
     o = int(statevars(ec_disotropic))
   end subroutine

   subroutine setmhdmagneticpressure(i)
     real, intent(in) :: i
     statevars(ec_muzero) = i
   end subroutine
   pure subroutine getmhdmagneticpressure(o)
     real, intent(out) :: o
     o = statevars(ec_muzero)
   end subroutine

   subroutine settfinish(i)
     real, intent(in) :: i
     statevars(ec_tfinish) = i
   end subroutine
   pure subroutine gettfinish(o)
     real, intent(out) :: o
     o = statevars(ec_tfinish)
   end subroutine

   subroutine setnpics(i)
     integer, intent(in) :: i
     if (i < 0) then
       call error('Number of snapshots cannot be negative', i, __FILE__, __LINE__)
     else
       statevars(ec_npics) = i
     end if
   end subroutine
   pure subroutine getnpics(o)
     integer, intent(out) :: o
     o = int(statevars(ec_npics))
   end subroutine

   subroutine sethfac(i)
     real, intent(in) :: i
     statevars(ec_hfac) = i
   end subroutine
   pure subroutine gethfac(o)
     real, intent(out) :: o
     o = statevars(ec_hfac)
   end subroutine

   subroutine setdtprint(i)
     real, intent(in) :: i
     statevars(ec_dtprint) = i
   end subroutine
   pure subroutine getdtprint(o)
     real, intent(out) :: o
     o = statevars(ec_dtprint)
   end subroutine

   subroutine setsilentmode(i)
     character (len=*), intent(in) :: i
     select case(i)
     case('yes')
       statevars(ec_silent) = 1
     case('no')
       statevars(ec_silent) = 0
     case('')
       statevars(ec_silent) = 0
     case default
       print *, 'There is no such silent mode: ', i
       stop
     end select
   end subroutine
   pure subroutine getsilentmode(o)
     integer, intent(out) :: o
     o = int(statevars(ec_silent))
   end subroutine

   subroutine setresolution(i)
     integer, intent(in) :: i
     statevars(ec_resolution) = i
   end subroutine
   pure subroutine getresolution(o)
     integer, intent(out) :: o
     o = int(statevars(ec_resolution))
   end subroutine

   subroutine setspacing(i)
     real, intent(in) :: i
     statevars(ec_spacing) = i
   end subroutine
   pure subroutine getspacing(o)
     real, intent(out) :: o
     o = statevars(ec_spacing)
   end subroutine

  subroutine setBorders(x1,x2,y1,y2,z1,z2,bs,pd)
  real, intent(in) :: &
   x1,x2,y1,y2,z1,z2,bs,pd
    statevars(ec_bordsize) = bs
    statevars(ec_padding) = pd
    statevars(ec_xmin) = x1
    statevars(ec_xmax) = x2
    statevars(ec_ymin) = y1
    statevars(ec_ymax) = y2
    statevars(ec_zmin) = z1
    statevars(ec_zmax) = z2
  end subroutine setBorders
  subroutine getBorders(x1,x2,y1,y2,z1,z2,bs,pd)
  real, intent(out) :: &
   x1,x2,y1,y2,z1,z2,bs,pd
    bs = statevars(ec_bordsize)
    pd = statevars(ec_padding)
    x1 = statevars(ec_xmin)
    x2 = statevars(ec_xmax)
    y1 = statevars(ec_ymin)
    y2 = statevars(ec_ymax)
    z1 = statevars(ec_zmin)
    z2 = statevars(ec_zmax)
  end subroutine

   subroutine setPartNumber(r, f)
     integer, optional, intent(in) :: r, f
     if (present(r)) statevars(ec_realpn) = r
     if (present(f)) statevars(ec_fixedpn) = f
   end subroutine
  subroutine getPartNumber(r, f)
    integer, optional, intent(out) :: r, f
    if (present(r)) r = int(statevars(ec_realpn))
    if (present(f)) f = int(statevars(ec_fixedpn))
  end subroutine

   subroutine setState(i)
     real, intent(in) :: i(ec_total)
     statevars(:) = i(:)
   end subroutine
   subroutine getState(o)
     real, intent(out) :: o(ec_total)
     o(:) = statevars(:)
   end subroutine

   subroutine setGamma(i)
     real, intent(in) :: i
     statevars(ec_gamma) = i
   end subroutine
   subroutine getGamma(o)
     real, intent(out) :: o
     o = statevars(ec_gamma)
   end subroutine

   subroutine setLastPrint(i)
     integer, intent(in) :: i
     statevars(ec_lastprint) = i
   end subroutine
   subroutine getLastPrint(o)
     integer, intent(out) :: o
     o = int(statevars(ec_lastprint))
   end subroutine

   subroutine setUseDumps(i)
     character (len=*), intent(in) :: i
     select case(i)
     case('yes')
       statevars(ec_usedumps) = 1
     case('no')
       statevars(ec_usedumps) = 0
     case('')
       statevars(ec_usedumps) = 1
     case default
       print *, 'There is no such mode for Use Dumps: ', i
       stop
     end select
   end subroutine
   pure subroutine getUseDumps(o)
     integer, intent(out) :: o
     o = int(statevars(ec_usedumps))
   end subroutine

   subroutine setProcess(it)
     character (len=*), intent(in) :: it
     select case(it)
     case('borderless')
       statevars(ec_process) = epc_borderless
     case('fullyperiodic')
       statevars(ec_process) = epc_fullyperiodic
     case('backcompatibility')
       statevars(ec_process) = epc_backcompatibility
     case default
       call error("The process is not set", it, __FILE__, __LINE__)
     end select
   end subroutine
   pure subroutine getProcess(ot)
     integer, intent(out) :: ot
     ot = int(statevars(ec_process))
   end subroutine

   subroutine splacement(itt, ott)
     character (len=*), intent(in) :: itt
     real, intent(out) :: ott

     select case(itt)
     case('uniform')
       ott = epl_uniform
     case('random')
       ott = epl_random
     case('closepacked')
       ott = epl_closepacked
     case default
       ott = e_none
       call error('There is no such initial particles distribution', itt, __FILE__, __LINE__)
     end select
   end subroutine

   subroutine printstate()
     use omp_lib
     use kernel_base,  only: kernelname

     implicit none

     real :: tmp
     character(len=*), parameter :: &
       bffs  = "(A,A50,ES13.4)"

     print *, '##############################################'
     print *, '#####'
     print*, "#   #"
     write(*,blockFormatInt) " #   #", "dim: ", int(statevars(ec_dim))

     if (int(statevars(ec_coordsys)) == 1) write(*,blockFormatStr) " #   #", "coordinate system: ", "Cartesian"
     if (int(statevars(ec_coordsys)) == 2) write(*,blockFormatStr) " #   #", "coordinate system: ", "Cylindric"

     select case(int(statevars(ec_eqs)))
     case(eeq_hydro)
       write(*,blockFormatStr) " #   #", "equations: ", "Hydro"
     case(eeq_diffusion)
       write(*,blockFormatStr) " #   #", "equations: ", "Diffusion"
     case(eeq_magnetohydro)
       write(*,blockFormatStr) " #   #", "equations: ", "Magneto Hydro"
     case(eeq_magnetohydrodiffusion)
       write(*,blockFormatStr) " #   #", "equations: ", "Magneto Hydro Diffusion"
     end select

     select case(int(statevars(ec_ics)))
     case(ett_sin3)
       write(*,blockFormatStr) " #   #", "initial conditions: ", "T = sin"
     case(ett_mti)
       write(*,blockFormatStr) " #   #", "initial conditions: ", "MTI Hopkins 2017"
     case(ett_shock12)
       write(*,blockFormatStr) " #   #", "initial conditions: ", "x < 0 => T = 1, x > 0 T = 2"
     case(ett_pulse)
       write(*,blockFormatStr) " #   #", "initial conditions: ", "T = gaussian pulse"
     case(ett_ring)
       write(*,blockFormatStr) " #   #", "initial conditions: ", "T = gaussian pulse in cylindring"
     case(ett_soundwave)
       write(*,blockFormatStr) " #   #", "initial conditions: ", "V,rho => soundwave"
     case(ett_hydroshock)
       write(*,blockFormatStr) " #   #", "initial conditions: ", "Hydroshock"
     case(ett_alfvenwave)
       write(*,blockFormatStr) " #   #", "initial conditions: ", "Alfven wave"
     case(ett_OTvortex)
       write(*,blockFormatStr) " #   #", "initial conditions: ", "Orszag-Tang vortex"
     end select

     if (int(statevars(ec_adden)) == 1) then
       write(*,blockFormatStr) " #   #", "advanced density: ", "yes"
     else
       write(*,blockFormatStr) " #   #", "advanced density: ", "no"
     end if

     if (int(statevars(ec_artts)) == 1) then
       write(*,blockFormatStr) " #   #", "artificail terms: ", "yes"
     else
       write(*,blockFormatStr) " #   #", "artificail terms: ", "no"
     end if

     select case(int(statevars(ec_ddw)))
     case(esd_fab)
       write(*,blockFormatStr) " #   #", "second derivative type: ", "Brookshaw"
     case(esd_fw)
       write(*,blockFormatStr) " #   #", "second derivative type: ", "New Brookshaw"
     case(esd_n2w)
       write(*,blockFormatStr) " #   #", "second derivative type: ", "Direct Derivative"
     case(esd_2nw_ds)
       write(*,blockFormatStr) " #   #", "second derivative type: ", "Two First Derivatives (-/+)"
     case(esd_2nw_sd)
       write(*,blockFormatStr) " #   #", "second derivative type: ", "Two First Derivatives (+/-)"
     end select
     write(*,blockFormatFltSci) " #   #", "print dt: ", statevars(ec_dtprint)
     write(*,blockFormatFltSci) " #   #", "stop time: ", statevars(ec_tfinish)
     write(*,blockFormatFlt) " #   #", "hfac: ", statevars(ec_hfac)

     if (int(statevars(ec_silent)) == 1) then
       write(*,blockFormatStr) " #   #", "silent mode: ", "yes"
     else
       write(*,blockFormatStr) " #   #", "silent mode: ", "no"
     end if
     if (int(statevars(ec_usedumps)) == 1) then
       write(*,blockFormatStr) " #   #", "use restore dumps: ", "yes"
     else
       write(*,blockFormatStr) " #   #", "use restore dumps: ", "no"
     end if
#ifdef _OPENMP
     write(*,blockFormatInt) " #   #", "OMP threads: ", int(omp_get_max_threads())
#endif
     if (int(statevars(ec_process)) == epc_fullyperiodic) then
       write(*,blockFormatStr) " #   #", "process: ", "fully periodic"
     else if (int(statevars(ec_process)) == epc_borderless) then
         write(*,blockFormatStr) " #   #", "process: ", "borderless"
     else
       write(*,blockFormatStr) " #   #", "process: ", "backcompatibility for hardcoded initial value"
     end if

     write(*,blockFormatInt) " #   #", "desired resolution: ", int(statevars(ec_resolution))

     write(*,blockFormatFlt2Sci) " #   #", "x in: [",statevars(ec_xmin)," :",statevars(ec_xmax)," ]"
     write(*,blockFormatFlt2Sci) " #   #", "y in: [",statevars(ec_ymin)," :",statevars(ec_ymax)," ]"
     write(*,blockFormatFlt2Sci) " #   #", "z in: [",statevars(ec_zmin)," :",statevars(ec_zmax)," ]"
     write(*,blockFormatFltSci) " #   #", "border size: ", statevars(ec_bordsize)

     write(*,blockFormatInt) " #   #", "real particles: ", int(statevars(ec_realpn))
     write(*,blockFormatInt) " #   #", "fixed particles: ", int(statevars(ec_fixedpn))
     write(*,blockFormatFlt) " #   #", "adiabatic gamma: ", statevars(ec_gamma)

     select case(int(statevars(ec_placement)))
     case(epl_uniform)
       write(*,blockFormatStr) " #   #", "initial particles placement: ", "Uniform"
     case(epl_random)
       write(*,blockFormatStr) " #   #", "initial particles placement: ", "Random"
     case(epl_closepacked)
       write(*,blockFormatStr) " #   #", "initial particles placement: ", "ClosePacked"
     end select

     if (statevars(ec_eqondiff)==1) then
       write(*,blockFormatStr) " #   #", "diffision equations: ", "yes"
     else
       write(*,blockFormatStr) " #   #", "diffision equations: ", "no"
     endif
     if (statevars(ec_eqonhydro)==1) then
       write(*,blockFormatStr) " #   #", "hydro equations: ", "yes"
     else
       write(*,blockFormatStr) " #   #", "hydro equations: ", "no"
     endif
     if (statevars(ec_eqonfluxlim)==1) then
       write(*,blockFormatStr) " #   #", "flux-limiter equations: ", "yes"
     else
       write(*,blockFormatStr) " #   #", "flux-limiter equations: ", "no"
     endif
     if (statevars(ec_eqonradexch)==1) then
       write(*,blockFormatStr) " #   #", "rad exchange equations: ", "yes"
     else
       write(*,blockFormatStr) " #   #", "rad exchange equations: ", "no"
     endif
     if (statevars(ec_eqonsts)==1) then
       if (statevars(ec_ststype)==0) then
         write(*,blockFormatStr) " #   #", "super time stepping: ", "auto"
       else
         write(*,blockFormatIntStr) " #   #", "super time stepping: ",&
          int(statevars(ec_stsfixeds))," stages"
       endif
     else
       write(*,blockFormatStr) " #   #", "super time stepping: ", "no"
     endif

     print*, "#   #"
     print*, "#####"
     print*, "##############################################"
   end subroutine printstate
end module
