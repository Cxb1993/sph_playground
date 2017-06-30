module n2ext
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase

  private
    ! real :: n2C(3) = (/ 0.123558082, 0.2994570731/pi, 0.236804709/pi /)
    ! sinc
    ! real :: n2C(3) = (/ 9.45261263832083, 7.29241786955192, 5.76692759942686 /)
    ! mgaus 031
    real :: n2C(3) = [0.0114824707771515, 0.00931249177952612, 0.00734547739346816]
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
      ! mgaus 031
      f  = -0.555723251459595*q**3 + 2.09531064900758*q**2 - 46.6863743125984*q  &
           - (-45.3101260774058*q + 41.4784899250259)*erf(1.7455062837353694*q  &
           - 1.5978981095815163)  &
           + 1.1397883441907*exp(q*(-3.04679218655966*q  &
           + 5.578282382086809)) + 40.2536188736974
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
      ! mgaus 031
      df  = -1.66716975437879*q**2 + 4.19062129801515*q + 45.3101260774058*erf(1.7455062837353694*q  &
            - 1.5978981095815163)  &
            - 46.6863743125984
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
      ! mgaus 031
      ddf  = -3.33433950875757*q + 4.19062129801515 + 89.2425038254715*exp(-3.04679218655966*(q  &
             - 0.91543532353377)**2)
    end if
  end subroutine
end module n2ext
