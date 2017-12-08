module args
  use errteylor,  only: setInfluenceCalc
  use state,      only: set_tasktype,&
                        setkerntype,&
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
    subroutine fillargs(dim, pspc1, resol, itype, ktype, dtype, errfname, dtout, npic, tfinish, sk, silent)
      real, intent(inout)                 :: pspc1, dtout, npic, tfinish, sk
      integer, intent(inout)              :: dim, silent, resol
      character (len=100), intent(inout)  :: itype, ktype, errfname, dtype
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
      ktype = ''
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

      print*, "# #"

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
          case('--kerneltype')
            ktype = adjustl(argval1)
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
      call setkerntype(ktype)
      dtout = tfinish / npic
      call set_difftype(dtype)

      print*, "# #"
      print *, "# #                dim:   ", dim
      if (coordsysstr /= '') then
        print *, "# #  coordinate system:   ", coordsysstr
      end if
      print *, "# #          task type:   ", itype
      if (initvart /= '') then
        call sinitvar(initvart)
        print *, "# #      init var type:   ", initvart
      end if
      call setAdvancedDensity(adden)
      print *, "# #   advanced density:   ", adden
      if (artts /= '') then
        print *, "# #   artificail terms:   ", artts
      end if
      print *, "# #           errfname:   ", errfname
      if (kerninflname /= '') then
        call setInfluenceCalc(kerninflname)
        print *, "# #     kern infl name:   ", kerninflname
      end if
      print *, "# #        kernel type:   ", ktype
      call getkernelname(kname)
      print *, "# #        kernel name:  ", kname
      write(*, "(A, F9.7)") " # #           print dt:   ", dtout
      write(*, "(A, F7.5)") " # #                  h:   ", sk
      write(*, "(A, A)") " # #           difftype:   ", dtype
      write(*, "(A, A)") " # #             silent:   ", silentstr
      if (pspc1 > 0) then
        write(*, "(A, F7.5)") " # #         desired dx:   x1=", pspc1
      end if
      if (resol > 0) then
        write(*, "(A, I5)") " # # desired resolution:   npx1=", resol
      end if

      call initkernel(dim)
    end subroutine fillargs
end module
