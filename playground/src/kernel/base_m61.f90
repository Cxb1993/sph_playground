module base_m61
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase, fwc

  private

    real :: n2C(3) = (/ 0.0250000245894439, 0.0419531421032309, 0.0716202677260738 /)
    real :: fwcl(3) = [4.*9., 4.*9.3125, 4.*9.6441]
    character (len=10) :: kernelname = ' M6/1 '
    real :: krad = 1.0, wCv, fwc
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

    if (q < 0.0) then
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    elseif (3*q < 1.0) then
      f  = 15.0*(-3.0*q + 1.0)**5 - 6.0*(-3.0*q + 2.0)**5 + (-3.0*q  &
+ 3.0)**5
    elseif (3*q < 2.0) then
      f  = -6.0*(-3.0*q + 2.0)**5 + (-3.0*q + 3.0)**5
    elseif (3*q < 3.0) then
      f  = (-3.0*q + 3.0)**5
    else
      f  = 0
    end if
  end subroutine

  pure subroutine kdf(q, df)
    real, intent(in)  :: q
    real, intent(out) :: df

    if (q < 0.0) then
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    elseif (3*q < 1.0) then
      df  = q*(-12150.0*q**3 + 9720.0*q**2 - 1080.0)
    elseif (3*q < 2.0) then
      df  = -15.0*(3.0*q - 3.0)**4 + 90.0*(3.0*q - 2.0)**4
    elseif (3*q < 3.0) then
      df  = -15.0*(3.0*q - 3.0)**4
    else
      df  = 0
    end if
  end subroutine

  pure subroutine kddf(q, ddf)
    real, intent(in)  :: q
    real, intent(out) :: ddf

    if (q < 0.0) then
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    elseif (3*q < 1.0) then
      ddf  = -48600.0*q**3 + 29160.0*q**2 - 1080.0
    else if (3*q < 2.0) then
      ddf  = 24300.0*q**3 - 43740.0*q**2 + 24300.0*q - 3780.0
    else if (3*q < 3.0) then
      ddf  = -180.0*(3.0*q - 3.0)**3
    else
      ddf  = 0
    end if
  end subroutine
end module
