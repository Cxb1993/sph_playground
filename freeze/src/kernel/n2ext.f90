module n2ext
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase

  private
    ! real :: n2C(3) = (/ 0.123558082, 0.2994570731/pi, 0.236804709/pi /)
    ! sinc
    ! real :: n2C(3) = (/ 9.45261263832083, 7.29241786955192, 5.76692759942686 /)
    ! unit
    real :: n2C(3) = [0.586533029501960, 0.436507677680717, 0.316192936870966]
    character (len=10) :: kernelname = ' genesis '
    real :: krad = 2.0, wCv
    integer :: dim
  contains

  subroutine setdimbase(d)
    integer, intent(in) :: d
    dim = d
    wCv = n2C(dim)
  end subroutine

  pure subroutine kf(q, f)
  ! subroutine kf(q, f)
    real, intent(in)  :: q
    real, intent(out) :: f

    if (q < 0.0) then
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative in f'
      end if
    else
      ! sinc
      ! f = 4*sin(pi*q/2)**4/(pi**5*q**4)
      ! unit
      if (q < 1.3333333333) then
        f = q*(0.5*q - 1.75) + 1.5
      else
        f = q*(-0.25*q**2 + 1.5*q - 3.0) + 2.0
      end if
    end if
  end subroutine

  pure subroutine kdf(q, df)
    real, intent(in)  :: q
    real, intent(out) :: df

    if (q < 0.0) then
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative in df'
      end if
    else
      ! sinc
      ! df = 8*(pi*q*cos(pi*q/2) - 2*sin(pi*q/2))*sin(pi*q/2)**3/(pi**5*q**5)
      ! unit
      if ( q < 1.333333333333 ) then
        df = q - 1.75
      else
        df = q*(-0.75*q + 3.0) - 3.0
      end if
    end if
  end subroutine

  pure subroutine kddf(q, ddf)
    real, intent(in)  :: q
    real, intent(out) :: ddf

    if (q < 0.0) then
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative in ddf'
      end if
    else
      ! sinc
      ! ddf  = 4*(2*pi**2*q**2*cos(pi*q) + pi**2*q**2 - 8*pi*q*sin(pi*q) - 10*cos(pi*q) + 10)*sin(pi*q/2)**2/(pi**5*q**6)
      ! unit
      if (q < 1.333333333333) then
        ddf = 1
      else
        ddf = -1.5*q + 3.0
      end if
    end if
  end subroutine
end module n2ext
