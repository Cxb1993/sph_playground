module diff
  use kernel,     only: nw
  use errprinter, only: error

  implicit none

  public :: c12nw, c2ddt, initdiffusion

  procedure(aic12nw), pointer :: c12nw => null()
  procedure(aic2ddt), pointer :: c2ddt => null()

  abstract interface
    subroutine aic12nw(r, ra, rb, dr, ha, hb, ma, fa, fb, dfdx)
      real, intent(out) :: dfdx(3)
      real, intent (in) :: r(3), ra(3), rb(3), dr, ha, hb, ma, fa, fb
    end subroutine aic12nw
  end interface

  abstract interface
    subroutine aic2ddt(r, ra, rb, dr, ha, hb, ma, fa, fb, dfdx)
      real, intent(out) :: dfdx(3)
      real, intent (in) :: r(3), ra(3), rb(3), dr, ha, hb, ma, fa, fb
    end subroutine aic2ddt
  end interface

contains
  subroutine initdiffusion()
    integer :: kt
    call getddwtype(kt)
    

    select case(kt)
    case(esd_n2w)
    case(esd_fab)
    case(esd_fw)
    case(esd_2nw_sd)
      c12nw => c12nwsum
    case(esd_2nw_ds)
      c12nw => c12nwdif
    case default
      call error("Wrong kernel type", kt, __FILE__, __LINE__)
    end select
  end subroutine initdiffusion

  pure subroutine c12nwdif(r, ra, rb, dr, ha, hb, ma, fa, fb, dfdx)
    real, intent(in) :: &
      r(3), ra(3), rb(3), dr, ha
    real, intent(out) :: &
      dtadx(3)
    real :: &
      nwa(3)

    call nw(r(:), ra(:), rb(:), dr, ha, nwa)
    dfdx(:) = mb*(fb - fa)*nwa(:)
  end subroutine c12nwdif

  pure subroutine c12nwsum(r, ra, rb, dr, ha, hb, ma, fa, fb, dfdx)
    real, intent(in) :: &
      r(3), ra(3), rb(3), dr, ha
    real, intent(out) :: &
      dtadx(3)
    real :: &
      nwa(3)

    call nw(r(:), ra(:), rb(:), dr, ha, nwa)
    call nw(r(:), ra(:), rb(:), dr, hb, nwb)
    dfdx(:) = mb*(fa*nwa(:) + fb*nmb(:))
  end subroutine c12nwsum

end module
