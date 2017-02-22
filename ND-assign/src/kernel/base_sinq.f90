module sinc
  use const

  implicit none

  public :: kf, kdf, kddf, knorm, krad, kernelname

  private

    real :: knorm(3) = (/ 0.123558082, 0.2994570731/pi, 0.236804709/pi /)
    real :: krad = 2.
    character (len=10) :: kernelname='sinc'

 contains

  subroutine kf(r, h, f)
    real, intent(in)  :: r, h
    real, intent(out) :: f
    real              :: q, term

    q = r / h
    if (q >= 2.) then
      f  = 0.
    else if (q >= 0.) then
      term = 0.5*pi*q
      f = (sin(term)**4)/(q**4)
    else
      print *, 'something went wrong, q =', q
      stop
    end if
  end subroutine kf

  subroutine kdf(r, h, df)
    real, intent(in)  :: r, h
    real, intent(out) :: df
    real              :: q, term

    q = r / h
    if (q >= 2.) then
      df = 0.
    else if (q >= 0.) then
      term = 0.5*pi*q
      df = 2./(q**4)*(sin(term)**3*cos(term)*pi - 2.*sin(term)**4/q)
    else
      print *, 'something went wrong, q =', q
      stop
    end if
  end subroutine kdf

  subroutine kddf(r, h, ddf)
    real, intent(in)  :: r, h
    real, intent(out) :: ddf
    real              :: q, term

    q = r / h
    if (q >= 2.) then
      ddf = 0.
    else if (q >= 0.) then
      term = 0.5*pi*q
      ddf = 1./(q**4)*(3.*sin(term)**2*cos(term)**2*pi**2 &
                  -16.*sin(term)**3*cos(term)*pi/q - sin(term)**4*pi**2 &
                  + 20.*(sin(term)**4)/q**2)
    else
      print *, 'something went wrong, q =', q
      stop
    end if
  end subroutine kddf
end module sinc
