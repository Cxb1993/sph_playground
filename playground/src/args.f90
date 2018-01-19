module args
  use const
  use errteylor,  only: setInfluenceCalc
  use state,      only: set_tasktype,&
                        setddwtype,&
                        set_difftype,&
                        sinitvar,&
                        setAdvancedDensity, &
                        setArtificialTerms, &
                        scoordsys
  use kernel,     only: initkernel, &
                        getkernelname

  implicit none

  public :: fillargs

  contains
    subroutine fillargs(dim, pspc1, resol, itype, ddwtype, dtype, errfname, dtout, npic, tfinish, sk, silent)
      real, intent(inout)                 :: pspc1, dtout, npic, tfinish, sk
      integer, intent(inout)              :: dim, silent, resol
      character (len=100), intent(inout)  :: itype, ddwtype, errfname, dtype
      character (len=100)                 :: kname

      integer                             :: numargs, curargnum
      character (len=100)                 :: argkey, argval1, silentstr, kerninflname, initvart,&
                                             adden, artts, coordsysstr
      real :: tmp

      dim   = 1
      itype = 'chi-laplace'
      pspc1 = 0
      resol = 0
      errfname = 'runresult.info'
      ddwtype = ''
      kname = ''
      tfinish = 1.
      npic = 200.
      sk = 1.2
      dtype = 'diff'
      silent = 0
      silentstr = 'no'
      kerninflname = ''
      initvart = ''
      adden = 'yes'
      artts = ''
      coordsysstr = ''

      print*, "#  #"

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
            itype = adjustl(argval1)
          case('--spacing')
            read(argval1, *) pspc1
          case('--resolution')
            read(argval1, *) tmp
            resol = int(tmp)
          case('--errfilename')
            errfname = adjustl(argval1)
          case('--kerninfluencefile')
            kerninflname = adjustl(argval1)
          case('--ddw')
            ddwtype = adjustl(argval1)
          case('--tfinish')
            read(argval1, *) tfinish
          case('--hfac')
            read(argval1, *) sk
          case('--difftype')
            dtype = adjustl(argval1)
          case('--silent')
            if (argval1 == "yes") then
              silent = 1
            else
              silent = 0
            end if
            silentstr = argval1
          case('--initvar')
            initvart = adjustl(argval1)
          case('--useadvanceddensity')
            adden = adjustl(argval1)
          case('--useartificialterms')
            artts = adjustl(argval1)
          case('--coordsys')
            coordsysstr = adjustl(argval1)
          case default
            print*, 'argument not found: ', argkey
            stop
          end select
          curargnum = curargnum + 1
        end do
      end if

      call scoordsys(coordsysstr)
      call setArtificialTerms(artts)
      call set_tasktype(itype)
      call setddwtype(ddwtype)
      dtout = tfinish / npic
      call set_difftype(dtype)

      print*, "#  #"
      write(*,blockFormatInt) " #  #", "dim:", dim
      if (coordsysstr /= '') then
        write(*,blockFormatStr) " #  #", "coordinate system: ", coordsysstr
      end if
      write(*,blockFormatStr) " #  #", "equations: ", itype
      if (initvart /= '') then
        call sinitvar(initvart)
        write(*,blockFormatStr) " #  #", "task type: ", initvart
      end if
      call setAdvancedDensity(adden)
      write(*,blockFormatStr) " #  #", "advanced density: ", adden
      if (artts /= '') then
        write(*,blockFormatStr) " #  #", "artificail terms: ", artts
      end if
      write(*,blockFormatStr) " #  #", "result file: ", errfname
      if (kerninflname /= '') then
        call setInfluenceCalc(kerninflname)
        write(*,blockFormatStr) " #  #", "kernel influence file name: ", kerninflname
      end if
      if (ddwtype /= '') then
        write(*,blockFormatStr) " #  #", "kernel type: ", ddwtype
      end if
      call getkernelname(kname)
      write(*,blockFormatStr) " #  #", "kernel name: ", kname
      write(*,blockFormatFlt) " #  #", "print dt: ", dtout
      write(*,blockFormatFlt) " #  #", "hfac: ", sk
      write(*,blockFormatStr) " #  #", "second derivative type: ", dtype
      write(*,blockFormatStr) " #  #", "silent mode: ", silentstr
      if (pspc1 > 0) then
        write(*,blockFormatFlt) " #  #", "desired spacing: ", pspc1
      end if
      if (resol > 0) then
        write(*,blockFormatInt) " #  #", "desired resolution:", resol
      end if

      call initkernel(dim)
    end subroutine fillargs
end module
