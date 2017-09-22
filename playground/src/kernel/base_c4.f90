module base_kernel
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase, fwc

  private

    real :: n2C(3)  = [ 3./4., 9./4./pi, 495./256./pi ]
    real :: fwcl(3) = [ 7.5, 9.428526, 9.75 ]
    character (len=20) :: kernelname = 'Wendland C4'
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
          f = (2. * q*q + 2.5 * q + 1.) * (1. - q/2.)**5
        else
          f = 0.
        end if
      else
        if (q < 2.) then
          f = (35./12. * q*q + 3. * q + 1.) * (1. - q/2.)**6
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
          df = (1. - 0.5*q)**5*(2.5 + 4.*q) - 2.5*(1. - 0.5*q)**4*(1. + 2.5*q + 2.*q**2)
        else
          df = 0.
        end if
      else
        if (q < 2.) then
          df = (1. - 0.5*q)**6*(3. + 5.833333333333333*q) - 3.*(1. - 0.5*q)**5*(1. + 3.*q + 2.9166666666666665*q**2)
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
          ddf = 4.*(1. - 0.5*q)**5 - 5.*(1. - 0.5*q)**4*(2.5 + 4.*q) + 5.*(1. - 0.5*q)**3*(1. + 2.5*q + 2.*q**2)
        else
          ddf = 0.
        end if
      else
        if (q < 2.) then
          ddf = 5.833333333333333*(1. - 0.5*q)**6 - 6.*(1. - 0.5*q)**5*(3. + &
                5.833333333333333*q) + 7.5*(1. - 0.5*q)**4*(1. + 3.*q + &
                2.9166666666666665*q**2)
        else
          ddf = 0.
        end if
      end if
    end if
  end subroutine
end module
