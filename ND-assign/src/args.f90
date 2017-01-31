module args
  use kernel, only:set_tasktype, set_kerntype, set_dim

  implicit none

  public :: fillargs

  private
    character (len=10) :: arg

  contains
    subroutine fillargs(dim, pspc1, pspc2, itype, ktype, errfname, dtout, npic, tfinish, sk)
      real, intent(inout)               :: pspc1, pspc2, dtout, npic, tfinish, sk
      integer, intent(inout)            :: dim
      character (len=40), intent(inout) :: itype, ktype, errfname

      call get_command_argument(1, arg)
      read(arg(:), fmt="(i5)") dim
      print *, "# #         dim:", dim
      call set_dim(dim)

      call get_command_argument(2, itype)
      print *, "# #   task type:   ", itype
      call set_tasktype(itype)

      pspc1 = 1.
      call get_command_argument(3, arg(:))
      read(arg(:), *) pspc1

      errfname = 'err'
      call get_command_argument(4, errfname)
      print *, "# #    errfname:   ", errfname

      ktype = 'fab'
      call get_command_argument(5, ktype)
      print *, "# #    ker.type:   ", ktype
      call set_kerntype(ktype)

      pspc2 = pspc1
      write(*, "(A, F7.5, A, F7.5)") " # #      set dx:   x1=", pspc1, "   x2=", pspc2
      call get_command_argument(6, arg(:))
      read(arg(:), *) tfinish
      npic = 200.
      dtout = tfinish / npic
      write(*, "(A, F9.7)") " # #    print dt:   ", dtout

      call get_command_argument(7, arg(:))
      read(arg(:), *) sk
      write(*, "(A, F7.5)") " # #           h:   ", sk
    end subroutine fillargs
end module args
