module args
  use errteylor,  only: setInfluenceCalc
  use state,      only: set_tasktype,&
                        set_kerntype,&
                        set_difftype,&
                        sinitvar,&
                        setAdvancedDensity, &
                        setArtificialTerms
  use kernel,     only: setdimkernel, &
                        getkernelname

  implicit none

  public :: fillargs

  contains
    subroutine fillargs(dim, pspc1, pspc2, itype, ktype, dtype, errfname, dtout, npic, tfinish, sk, silent)
      real, intent(inout)                 :: pspc1, pspc2, dtout, npic, tfinish, sk
      integer, intent(inout)              :: dim, silent
      character (len=100), intent(inout)  :: itype, ktype, errfname, dtype
      character (len=100)                 :: kname

      integer                             :: numargs, curargnum
      character (len=100)                 :: argkey, argval1, silentstr, kerninflname, initvart,&
                                             adden, artts

      dim   = 1
      itype = 'chi-laplace'
      pspc1 = 1.
      errfname = 'runresult.info'
      ktype = 'fab'
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
      artts = 'yes'


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
          case('--tasktype')
            itype = adjustl(argval1)
          case('--spacing')
            read(argval1, *) pspc1
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
          case default
            print*, 'argument not found: ', argkey
            stop
          end select
          curargnum = curargnum + 1
        end do
      end if

      call setdimkernel(dim)
      call set_tasktype(itype)
      pspc2 = pspc1
      call set_kerntype(ktype)
      dtout = tfinish / npic
      call set_difftype(dtype)

      print *, "# #              dim:", dim
      print *, "# #        task type:   ", itype
      if (initvart /= '') then
        call sinitvar(initvart)
        print *, "# #    init var type:   ", initvart
      end if
      call setAdvancedDensity(adden)
      print *, "# # advanced density:   ", adden
      call setArtificialTerms(artts)
      print *, "# # artificail terms:   ", artts
      print *, "# #         errfname:   ", errfname
      if (kerninflname /= '') then
        call setInfluenceCalc(kerninflname)
        print *, "# #   kern infl name:   ", kerninflname
      end if
      print *, "# #         ker.type:   ", ktype
      call getkernelname(kname)
      print *, "# #         ker.name:   ", kname
      write(*, "(A, F9.7)") " # #         print dt:   ", dtout
      write(*, "(A, F7.5)") " # #                h:   ", sk
      write(*, "(A, A)") " # #         difftype:   ", dtype
      write(*, "(A, A)") " # #           silent:   ", silentstr
      write(*, "(A, F7.5, A, F7.5)") " # #           set dx:   x1=", pspc1, "   x2=", pspc2
    end subroutine fillargs
end module
