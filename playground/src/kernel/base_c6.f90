module base_kernel
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase, fwc

  private

    real :: n2C(3)  = [ 55./64., 78./28./pi, 1365./512./pi ]
    real :: fwcl(3) = [ 9.75, 11.66665, 12. ]
    character (len=20) :: kernelname = 'Wendland C6'
    real :: krad = 2.0, wCv, fwc
    integer :: dim
  contains

  subroutine setdimbase(d)
    integer, intent(in) :: d
    dim = d
    wCv = n2C(dim)
    fwc = fwcl(dim)
  end subroutine

  pure subroutine kf(q, f)
    real, intent(in)  :: q
    real, intent(out) :: f

    if (q < 0.) then
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    else
      if (dim == 1) then
        if (q < 2.) then
          f = (1. - q/2.)**7 * (21./8. * q*q*q + 19./4. * q*q + 3.5 * q + 1.)
        else
          f = 0.
        end if
      else
        if (q < 2.) then
          f = (1. - q/2.)**8 * (4. * q*q*q + 6.25 * q*q + 4. * q + 1.)
        else
          f = 0.
        end if
      end if
    end if
  end subroutine

  pure subroutine kdf(q, df)
    real, intent(in)  :: q
    real, intent(out) :: df

    if (q < 0.) then
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    else
      if (dim == 1) then
        if (q < 2.) then
          df = (1. - 0.5*q)**7*(3.5 + 9.5*q + 7.875*q**2) - 3.5*(1. - 0.5*q)**6*(1. + 3.5*q + 4.75*q**2 + 2.625*q**3)
        else
          df = 0.
        end if
      else
        if (q < 2.) then
          df = (1. - 0.5*q)**8*(4. + 12.5*q + 12.*q**2) - 4.*(1. - 0.5*q)**7*(1. + 4.*q + 6.25*q**2 + 4.*q**3)
        else
          df = 0.
        end if
      end if
    end if
  end subroutine

  pure subroutine kddf(q, ddf)
    real, intent(in)  :: q
    real, intent(out) :: ddf

    if (q < 0.) then
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    else
      if (dim == 1) then
        if (q < 2.) then
          ddf = (1. - 0.5*q)**7*(9.5 + 15.75*q) - 7.*(1. - 0.5*q)**6*(3.5 + 9.5*q + 7.875*q**2) + &
                10.5*(1. - 0.5*q)**5*(1. + 3.5*q + 4.75*q**2 + 2.625*q**3)
        else
          ddf = 0.
        end if
      else
        if (q < 2.) then
          ddf = (1. - 0.5*q)**8*(12.5 + 24.*q) - 8.*(1. - 0.5*q)**7*(4. + 12.5*q + &
                12.*q**2) + 14.*(1. - 0.5*q)**6*(1. + 4.*q + 6.25*q**2 + 4.*q**3)
        else
          ddf = 0.
        end if
      end if
    end if
  end subroutine
end module
