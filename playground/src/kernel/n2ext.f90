module n2ext
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase

  private

    real :: n2C(3) = (/ 6.55410349636388, 5.68497842356648, 4.73825676824363 /)
    character (len=10) :: kernelname = ' m4-- '
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
    elseif (q < 1.0) then
      f  = 0.0833335*q**2 - 0.22777878889*q + 0.159524
    elseif (q < 2.0) then
      f  = -0.00396825*q**7 + 0.0333333*q**6 - 0.1*q**5 + 0.111111*q**4  &
- 0.177778*q + 0.152381
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
    elseif (q < 1.0) then
      df  = 0.166667*q - 0.22777878889
    elseif (q < 2.0) then
      df  = -0.02777775*q**6 + 0.1999998*q**5 - 0.5*q**4 + 0.444444*q**3  &
- 0.177778
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
    elseif (q < 1.0) then
      ddf  = 0.1666670000000
    else if (q < 2.0) then
      ddf  = q**2*(-0.1666665*q**3 + 0.999999*q**2 - 2.0*q + 1.333332)
    else
      ddf  = 0
    end if
  end subroutine
end module n2ext

