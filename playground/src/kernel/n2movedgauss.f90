module n2movedgauss
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase

  private

    real :: n2C(3) = (/ 0.0114685539154613, 0.00930942111372185, 0.00735105715323899 /)
    character (len=10) :: kernelname = ' moved gauss '
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
    elseif (q < 2.0) then
      f  = -0.55976281236227*q**3 + 2.16373843962652*q**2 - 46.8240032936104*q  &
- (-45.2069218327059*q + 41.4282147429831)*erf(1.756941181815487*q  &
- 1.6100838902591093)  &
+ 1.08647749179408*exp(q*(-3.0868423163592*q  &
+ 5.657645385947832)) + 40.4459781700414
    else
      f = .0
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
    elseif (q < 2.0) then
      df  = -1.67928843708681*q**2 + 4.32747687925304*q + 45.2069218327059*erf(1.756941181815487*q  &
- 1.6100838902591093)  &
- 46.8240032936104
    else
      df = .0
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
    elseif (q < 2.0) then
      ddf  = -3.35857687417362*q + 4.32747687925304 + 89.6225339017062*exp(-3.0868423163592*(q  &
- 0.916413085949396)**2)
    else
      ddf = .0
    end if
  end subroutine
end module
