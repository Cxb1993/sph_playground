module state

  use const
  use errprinter, only: error

  implicit none

  public :: set_equations, get_equations, &
            setddwtype, getddwtype, getdim, setdim, sinitvar, &
            ginitvar, setAdvancedDensity, getAdvancedDensity,&
            setArtificialTerms, getArtificialTerms, &
            scoordsys, gcoordsys, &
            setdiffconductivity, getdiffconductivity, &
            setdiffisotropic, getdiffisotropic, &
            setmhdmagneticpressure, getmhdmagneticpressure, &
            printstate, setresultfile, getresultfile,&
            setkerninflfilename, getkerninflfilename,&
            setnpics, settfinish, gethfac, gettfinish,&
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
            setArtTermCond, getArtTermCond

  private
  save
    real :: statevars(ec_total)
    character (len=100) :: resfilename, kerninflname

  contains
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
        eqSet(eqs_hydro) = 1
        eqSet(eqs_magneto) = 1
      case(eeq_diffusion)
        eqSet(eqs_diff) = 1
      case(eeq_hyrad)
        eqSet(eqs_fld) = 1
        eqSet(eqs_hydro) = 1
      case(eeq_magnetohydrodiffusion)
        eqSet(eqs_hydro) = 1
        eqSet(eqs_magneto) = 1
        eqSet(eqs_diff) = 1
      case(eeq_hydrodiffusion)
        eqSet(eqs_hydro) = 1
        eqSet(eqs_diff) = 1
      case default
        call error("Equation is not set", int(statevars(ec_eqs)), __FILE__, __LINE__)
      end select
    end subroutine getEqComponent

    subroutine setddwtype(itt)
      character (len=*), intent(in) :: itt
      select case(itt)
      case('n2w')
        statevars(ec_ddw) = esd_n2w
      case('fab')
        statevars(ec_ddw) = esd_fab
      case('2nw_ds')
        statevars(ec_ddw) = esd_2nw_ds
      case('2nw_sd')
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

   subroutine sinitvar(itt)
     character (len=*), intent(in) :: itt
     select case(itt)
     case('sin3')
       statevars(ec_ics) = ett_sin3
     case('mti')
       statevars(ec_ics) = ett_mti
     case('mtilowres')
       statevars(ec_ics) = ett_mtilowres
     case('shock12')
       statevars(ec_ics) = ett_shock12
     case('pulse')
       statevars(ec_ics) = ett_pulse
     case('ring')
       statevars(ec_ics) = ett_ring
     case('soundwave')
       statevars(ec_ics) = ett_soundwave
     case('hydroshock')
       statevars(ec_ics) = ett_hydroshock
     case('alfvenwave')
       statevars(ec_ics) = ett_alfvenwave
     case('otvortex')
       statevars(ec_ics) = ett_OTvortex
     case('boilingtank')
       statevars(ec_ics) = ett_boilingtank
     case('fld_gauss')
       statevars(ec_ics) = ett_fld_gauss
     case default
       statevars(ec_ics) = e_none
       call error('There is no such initial variable', itt, __FILE__, __LINE__)
     end select
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

   subroutine setresultfile(i)
     character (len=*), intent(in) :: i
     if (i /= '') then
       resfilename = i
     else
       resfilename = 'result.info'
     end if
   end subroutine
   pure subroutine getresultfile(o)
     character (len=*), intent(out) :: o
     o = trim(resfilename)
   end subroutine

   subroutine setkerninflfilename(i)
     character (len=*), intent(in) :: i
     if (i /= '') then
       kerninflname = i
     else
       kerninflname = 'influence.info'
     end if
   end subroutine
   pure subroutine getkerninflfilename(o)
     character (len=*), intent(out) :: o
     o = trim(kerninflname)
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
   end subroutine
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
   ! pure subroutine getPartNumber(r, f)
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
     case('relaxation')
       statevars(ec_process) = epc_relaxation
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


   subroutine printstate()
     use kernel_base,  only: kernelname

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

     write(*,blockFormatStr) " #   #", "result file: ", trim(resfilename)
     write(*,blockFormatStr) " #   #", "kernel influence file name: ", trim(kerninflname)

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

     write(*,blockFormatStr) " #   #", "kernel name: ", kernelname
     write(*,blockFormatFlt) " #   #", "print dt: ", statevars(ec_dtprint)
     write(*,blockFormatFlt) " #   #", "stop time: ", statevars(ec_tfinish)
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
     if (int(statevars(ec_process)) == epc_relaxation) then
       write(*,blockFormatStr) " #   #", "process: ", "relaxation"
     else
       write(*,blockFormatStr) " #   #", "process: ", "backcompatibility"
     end if

write(*,blockFormatInt) " #   #", "desired resolution:", int(statevars(ec_resolution))

     print *, '#   #                            x in [',statevars(ec_xmin),":",statevars(ec_xmax),"]"
     print *, '#   #                            y in [',statevars(ec_ymin),":",statevars(ec_ymax),"]"
     print *, '#   #                            z in [',statevars(ec_zmin),":",statevars(ec_zmax),"]"
     write(*,blockFormatFlt) " #   #", "border size: ", statevars(ec_bordsize)

     write(*,blockFormatInt) " #   #", "real particles: ", int(statevars(ec_realpn))
     write(*,blockFormatInt) " #   #", "fixed particles: ", int(statevars(ec_fixedpn))
     write(*,blockFormatFlt) " #   #", "adiabatic gamma: ", statevars(ec_gamma)


     print*, "#   #"
     print*, "#####"
     print*, "##############################################"
   end subroutine printstate
end module
