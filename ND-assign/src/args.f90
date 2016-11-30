module args
  implicit none

  public :: fillargs

  private
    character (len=10) :: arg

  contains
    subroutine fillargs(dim, pspc1, pspc2, itype, ktype, errfname, dtout, npic, tfinish)
      real, intent(inout)               :: pspc1, pspc2, dtout, npic, tfinish
      integer, intent(inout)            :: dim
      character (len=40), intent(inout) :: itype, ktype, errfname

      call get_command_argument(1, arg)
      read(arg(:), fmt="(i5)") dim
      print *, "#         dim:", dim

      call get_command_argument(2, itype)
      print *, "#   task type:   ", itype

      pspc1 = 1.
      call get_command_argument(3, arg(:))
      read(arg(:), *) pspc1

      errfname = 'err'
      call get_command_argument(4, errfname)
      print *, "#    errfname:   ", errfname

      ktype = 'fab'
      call get_command_argument(5, ktype)
      print *, "#    ker.type:   ", ktype

      pspc2 = pspc1
      write(*, "(A, F7.5, A, F7.5)") " #      set dx:   x1=", pspc1, "   x2=", pspc2
      call get_command_argument(6, arg(:))
      read(arg(:), *) tfinish

      npic = 200.
      dtout = tfinish / npic
      write(*, "(A, F9.7)") " #    print dt:   ", dtout
    end subroutine fillargs
end module args
