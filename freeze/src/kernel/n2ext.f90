module n2ext
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase

  private

    real :: n2C(3) = (/ 0.999994577452564, 0.682084687280554, 0.477457160036309 /)
    character (len=10) :: kernelname = ' zipM6 '
    real :: krad = 2.0, wCv
    integer :: dim
  contains

  subroutine setdimbase(d)
    integer, intent(in) :: d
    dim = d
    wCv = n2C(dim)
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
    elseif (q <= 1.0) then
      f  = 0.5*q**3 - 1.0*q**2 + 0.6666666666667
    elseif (q <= 2.0) then
      f  = -0.1666666666667*(q - 2.0)**3
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
    elseif (q <= 1.0) then
      df  = q*(1.5*q - 2.0)
    elseif (q <= 2.0) then
      df  = -0.5*(q - 2.0)**2
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
    elseif (q <= 1.0) then
      ddf  = 3.0*q - 2.0
    else if (q <= 2.0) then
      ddf  = -1.0*q + 2.0
    else
      ddf  = 0
    end if
  end subroutine
end module
