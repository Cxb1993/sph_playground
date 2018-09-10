module args
  use const
  use state,      only: set_equations,&
                        setddwtype,&
                        sinitvar,&
                        setAdvancedDensity, &
                        setArtificialTerms, &
                        scoordsys,&
                        setresultfile,&
                        setkerninflfilename,&
                        setnpics,&
                        settfinish,&
                        setresolution,&
                        setspacing,&
                        sethfac,&
                        setsilentmode,&
                        setLastPrint,&
                        setUseDumps,&
                        setdim,&
                        setProcess
use errprinter,   only: warning
  implicit none

  public :: fillargs

  contains
    subroutine fillargs()
      real :: &
        pspc1, tfinish, hfac, dimf, npicf
      integer :: &
        numargs, curargnum, npic, dim, silent, resol
      character (len=100) :: &
        eqs, resultfile, ddwtype,&
        argkey, argval1, silentstr, kerninflname, initvart,&
        adden, artts, coordsysstr, usedumps, process
      real :: tmp

      dim = 0
      dimf = 0.0
      eqs = ''
      pspc1 = 0
      resol = 0
      resultfile = ''
      ddwtype = ''
      tfinish = 1.
      npic = 0
      npic = 0.0
      hfac = 1.2
      silent = 0
      silentstr = ''
      kerninflname = ''
      initvart = ''
      adden = ''
      artts = ''
      coordsysstr = ''
      usedumps = ''
      process = ''

      numargs = command_argument_count()
      if ( numargs > 0 )then
        curargnum = 1
        do while (curargnum < numargs)
          call get_command_argument(curargnum, argkey)
          curargnum = curargnum + 1
          call get_command_argument(curargnum, argval1)
          select case(adjustl(argkey))
          case('--dim')
            read(argval1, *) dimf
            dim = int(dimf)
          case('--equations')
            eqs = adjustl(argval1)
          case('--spacing')
            read(argval1, *) pspc1
          case('--nsnapshots')
            read(argval1, *) npicf
            npic = int(npicf)
          case('--resolution')
            read(argval1, *) tmp
            resol = int(tmp)
          case('--resultfile')
            resultfile = adjustl(argval1)
          case('--influencefile')
            kerninflname = adjustl(argval1)
          case('--ddw')
            ddwtype = adjustl(argval1)
          case('--tfinish')
            read(argval1, *) tfinish
          case('--hfac')
            read(argval1, *) hfac
          case('--silent')
            silentstr = adjustl(argval1)
          case('--ics')
            initvart = adjustl(argval1)
          case('--adden')
            adden = adjustl(argval1)
          case('--artts')
            artts = adjustl(argval1)
          case('--coordsys')
            coordsysstr = adjustl(argval1)
          case('--usedumps')
            usedumps = adjustl(argval1)
          case('--process')
            process = adjustl(argval1)
          case default
            call warning("Argument not found", argkey, __FILE__, __LINE__)
          end select
          curargnum = curargnum + 1
        end do
      end if

      if (npic /= 0) then
        call setnpics(npic)
      end if
      if (tfinish /= 0.0) then
        call settfinish(tfinish)
      end if
      if (resultfile /= '') then
        call setresultfile(resultfile)
      end if
      if (kerninflname /= '') then
        call setkerninflfilename(kerninflname)
      end if
      if (initvart /= '') then
        call sinitvar(initvart)
      end if
      if (coordsysstr /= '') then
        call scoordsys(coordsysstr)
      end if
      if (artts /= '') then
        call setArtificialTerms(artts)
      end if
      if (eqs /= '') then
        call set_equations(eqs)
      end if
      if (ddwtype /= '') then
        call setddwtype(ddwtype)
      end if
      if (hfac /= 1.2) then
        call sethfac(hfac)
      end if
      if (silentstr /= '') then
        call setsilentmode(silentstr)
      end if
      if (pspc1 /= 0) then
        call setspacing(pspc1)
      end if
      if (resol /= 0) then
        call setresolution(resol)
      end if
      if (adden /= '') then
        call setAdvancedDensity(adden)
      end if
      if (usedumps /= '') then
        call setUseDumps(usedumps)
      end if
      if (dim /= 0) then
        call setdim(dim)
      end if
      if (process /= '') then
        call setProcess(process)
      end if
      call setLastPrint(0)
    end subroutine fillargs
end module
