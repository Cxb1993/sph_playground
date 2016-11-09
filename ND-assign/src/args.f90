module args
  implicit none

  public :: fillargs

  private
    character (len=10) :: arg

  contains
    subroutine fillargs(dim, pspc1, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc2,&
                        itype, ktype, errfname, dtout, npic, tfinish)
      real, intent(inout)               :: pspc1, pspc2, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2,&
                                           dtout, npic, tfinish
      integer, intent(inout)            :: dim
      character (len=40), intent(inout) :: itype, ktype, errfname

      call get_command_argument(1, arg)
      read(arg(:), fmt="(i5)") dim
      print *, "#       dim:", dim

      call get_command_argument(2, itype)
      print *, "# task type:   ", itype

      pspc1 = 1.
      call get_command_argument(3, arg(:))
      read(arg(:), *) pspc1

      errfname = 'err'
      call get_command_argument(4, errfname)
      print *, "#  errfname:   ", errfname

      ktype = 'fab'
      call get_command_argument(5, ktype)
      print *, "#  ker.type:   ", ktype

      pspc2 = pspc1
      write(*, "(A, F7.5, A, F7.5)") " # p.spacing:   x1=", pspc1, "   x2=", pspc2
      brdx1 = -10.
      brdx2 = 10.
      if (dim.gt.1) then
        brdy1 = - pspc1 * 4.5
        brdy2 = pspc2 * 4.5
      else
        brdy1 = 0.
        brdy2 = 0.
      end if
      if (dim.eq.3) then
        brdz1 = - pspc1 * 4.5
        brdz2 = pspc2 * 4.5
      else
        brdz1 = 0.
        brdz2 = 0.
      end if

      tfinish = 5.
      npic = 200.
      dtout = tfinish / npic
      write(*, "(A, F9.7)") " #  print dt:   ", dtout
    end subroutine fillargs

end module args
