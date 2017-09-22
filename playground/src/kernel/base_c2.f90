module base_kernel
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase, fwc

  private

    real :: n2C(3)  = [ 5./8., 7./4./pi, 21./16./pi ]
    real :: fwcl(3) = [ 5.25, 7.2, 7.5 ]
    character (len=20) :: kernelname = ' Wendland C2 '
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
          f = (1.5 * q + 1.) * (1. - q/2.)**3
        else
          f = 0.
        end if
      else
        if (q < 2.) then
          f = (2. * q + 1.) * (1. - q/2.)**4
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
          df = 1.5*(1. - 0.5*q)**3 - 1.5*(1. - 0.5*q)**2*(1. + 1.5*q)
        else
          df = 0.
        end if
      else
        if (q < 2.) then
          df = 2.*(1. - 0.5*q)**4 - 2.*(1. - 0.5*q)**3*(1. + 2.*q)
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
          ddf = -4.5*(1. - 0.5*q)**2 + 1.5*(1. - 0.5*q)*(1. + 1.5*q)
        else
          ddf = 0.
        end if
      else
        if (q < 2.) then
          ddf = -8.*(1. - 0.5*q)**3 + 3.*(1. - 0.5*q)**2*(1. + 2.*q)
        else
          ddf = 0.
        end if
      end if
    end if
  end subroutine
end module
