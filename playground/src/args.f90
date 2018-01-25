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
                        setUseDumps
  use kernel,     only: initkernel, &
                        getkernelname

  implicit none

  public :: fillargs

  contains
    subroutine fillargs()
      real :: &
        pspc1, dtout, tfinish, hfac
      integer :: &
        numargs, curargnum, npic,dim, silent, resol
      character (len=100) :: &
        eqs, resultfile, ddwtype,&
        argkey, argval1, silentstr, kerninflname, initvart,&
        adden, artts, coordsysstr, usedumps
      real :: tmp

      dim   = 1
      eqs = ''
      pspc1 = 0
      resol = 0
      resultfile = ''
      ddwtype = ''
      tfinish = 1.
      npic = 200
      hfac = 1.2
      silent = 0
      silentstr = ''
      kerninflname = ''
      initvart = ''
      adden = ''
      artts = ''
      coordsysstr = ''
      usedumps = ''

      numargs = command_argument_count()
      if ( numargs > 0 )then
        curargnum = 1
        do while (curargnum < numargs)
          call get_command_argument(curargnum, argkey)
          curargnum = curargnum + 1
          call get_command_argument(curargnum, argval1)
          select case(adjustl(argkey))
          case('--dim')
            read(argval1, fmt="(i5)") dim
          case('--equations')
            eqs = adjustl(argval1)
          case('--spacing')
            read(argval1, *) pspc1
          case('--resolution')
            read(argval1, *) tmp
            resol = int(tmp)
          case('--resultfile')
            resultfile = adjustl(argval1)
          case('--kerninfluencefile')
            kerninflname = adjustl(argval1)
          case('--ddw')
            ddwtype = adjustl(argval1)
          case('--tfinish')
            read(argval1, *) tfinish
          case('--hfac')
            read(argval1, *) hfac
          case('--silent')
            silentstr = adjustl(argval1)
          case('--initvar')
            initvart = adjustl(argval1)
          case('--useadvanceddensity')
            adden = adjustl(argval1)
          case('--useartificialterms')
            artts = adjustl(argval1)
          case('--coordsys')
            coordsysstr = adjustl(argval1)
          case('--usedumps')
            usedumps = adjustl(argval1)
          case default
            print*, 'argument not found: ', argkey
            stop
          end select
          curargnum = curargnum + 1
        end do
      end if

      call setnpics(npic)
      call settfinish(tfinish)
      call setresultfile(resultfile)
      call setkerninflfilename(kerninflname)
      call sinitvar(initvart)
      call scoordsys(coordsysstr)
      call setArtificialTerms(artts)
      call set_equations(eqs)
      call setddwtype(ddwtype)
      dtout = tfinish / npic
      call initkernel(dim)
      call sethfac(hfac)
      call setsilentmode(silentstr)
      call setspacing(pspc1)
      call setresolution(resol)
      call setAdvancedDensity(adden)
      call setLastPrint(0)
      call setUseDumps(usedumps)
    end subroutine fillargs
end module
