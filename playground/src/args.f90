module args
  use errteylor,  only: setInfluenceCalc
  use kernel,     only: set_tasktype,&
                        set_kerntype,&
                        set_dim,&
                        set_difftype

  implicit none

  public :: fillargs

  contains
    ! ARGS Dim, TaskType, Spacing, ErrFileName, KernelType, tFinish, hfac
    subroutine fillargs(dim, pspc1, pspc2, itype, ktype, dtype, errfname, dtout, npic, tfinish, sk, silent)
      real, intent(inout)               :: pspc1, pspc2, dtout, npic, tfinish, sk
      integer, intent(inout)            :: dim, silent
      character (len=40), intent(inout) :: itype, ktype, errfname, dtype

      integer                           :: numargs, curargnum
      character (len=40)                :: argkey, argval, silentstr, kerninflname

      dim   = 1
      itype = 'chi-laplace'
      pspc1 = 1.
      errfname = 'runresult.info'
      ktype = 'fab'
      tfinish = 1.
      npic = 200.
      sk = 1.2
      dtype = 'diff'
      silent = 0
      silentstr = 'no'
      kerninflname = ''

      numargs = command_argument_count()
      if ( numargs > 0 )then
        curargnum = 1
        do while (curargnum < numargs)
          call get_command_argument(curargnum, argkey)
          curargnum = curargnum + 1
          call get_command_argument(curargnum, argval)
          select case(adjustl(argkey))
          case('--dim')
            read(argval, fmt="(i5)") dim
          case('--tasktype')
            itype = adjustl(argval)
          case('--spacing')
            read(argval, *) pspc1
          case('--errfilename')
            errfname = adjustl(argval)
          case('--kerninfluencefile')
            kerninflname = adjustl(argval)
          case('--kerneltype')
            ktype = adjustl(argval)
          case('--tfinish')
            read(argval, *) tfinish
          case('--hfac')
            read(argval, *) sk
          case('--difftype')
            dtype = adjustl(argval)
          case('--silent')
            if (argval == "yes") then
              silent = 1
            else
              silent = 0
            end if
            silentstr = argval
          case default
            print*, 'argument not found: ', argkey
            stop
          end select
          ! print*, argkey, argval
          curargnum = curargnum + 1
        end do
      end if

      call set_dim(dim)
      call set_tasktype(itype)
      pspc2 = pspc1
      call set_kerntype(ktype)
      dtout = tfinish / npic
      call set_difftype(dtype)

      print *, "# #            dim:", dim
      print *, "# #      task type:   ", itype
      print *, "# #       errfname:   ", errfname
      if (kerninflname /= '') then
        call setInfluenceCalc(kerninflname)
        print *, "# # kern infl name:   ", kerninflname
      end if
      print *, "# #       ker.type:   ", ktype
      write(*, "(A, F9.7)") " # #       print dt:   ", dtout
      write(*, "(A, F7.5)") " # #              h:   ", sk
      write(*, "(A, A)") " # #       difftype:   ", dtype
      write(*, "(A, A)") " # #         silent:   ", silentstr
      write(*, "(A, F7.5, A, F7.5)") " # #         set dx:   x1=", pspc1, "   x2=", pspc2
    end subroutine fillargs
end module args
