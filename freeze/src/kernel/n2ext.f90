module n2ext
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase

  private
    ! real :: n2C(3) = (/ 0.123558082, 0.2994570731/pi, 0.236804709/pi /)
    real :: n2C(3) = (/ 9.45261263832083, 7.29241786955192, 5.76692759942686 /)
    character (len=10) :: kernelname = ' genesis '
    real :: krad = 2.0, wCv
    integer :: dim
  contains

  subroutine setdimbase(d)
    integer, intent(in) :: d
    dim = d
    wCv = n2C(dim)
  end subroutine

  pure subroutine kf(r, h, f)
    real, intent(in)  :: r, h
    real, intent(out) :: f
    real              :: q

    q = r / h
    if (q < 0.0) then
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    else
      f = 4*sin(pi*q/2)**4/(pi**5*q**4)
    end if
  end subroutine

  pure subroutine kdf(r, h, df)
    real, intent(in)  :: r, h
    real, intent(out) :: df
    real              :: q

    q = r / h
    if (q < 0.0) then
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    else
      df = 8*(pi*q*cos(pi*q/2) - 2*sin(pi*q/2))*sin(pi*q/2)**3/(pi**5*q**5)
    end if
    df = df / q
  end subroutine

  pure subroutine kddf(r, h, ddf)
    real, intent(in)  :: r, h
    real, intent(out) :: ddf
    real              :: q

    q = r / h
    if (q < 0.0) then
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    else
      ddf  = 4*(2*pi**2*q**2*cos(pi*q) + pi**2*q**2 - 8*pi*q*sin(pi*q) - 10*cos(pi*q) + 10)*sin(pi*q/2)**2/(pi**5*q**6)
    end if
  end subroutine
end module n2ext
