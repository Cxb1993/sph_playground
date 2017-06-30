module n2q2m4
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase

  private

    real :: n2C(3) = (/ 6.61200174725079, 5.58719034631753, 4.60567114916346 /)
    character (len=40) :: kernelname = ' Int(Int(q^{2.11}*M4)) '
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
      f  = -0.1630499003148*q + 0.001569029064886*q**4.11*(3.0*q**3  &
- 8.0*q**2 + 8.0) + 0.0005517752211516*q**4.11*(3.0*q**3  &
- 6.0*q**2 + 4.0) + 0.03739693721433*q**4.11  &
- 0.01616571128374*q**6.11 + 0.005147159653345*q**7.11  &
+ 0.1470241578171
    elseif (q <= 2.0) then
      f  = -0.1217108669236*q - 0.0005517752211516*q**4.11*(q  &
- 2.0)**3 - 5.873741473634e-5*q**4.11*(q - 2.0)**2*(26.7126*q  &
- 106.8504) + 0.02819395907344*q**4.11  &
- 0.01409697953672*q**5.11 + 0.00352424488418*q**6.11  &
- 0.000352424488418*q**7.11 - 5.873741473634e-5*(2.11*q  &
- 4.22)*(188.0*q**4.11 - 80.0*q**5.11 +  &
11.0*q**6.11 + 414.4236428177) + 0.04868434676946
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
      df  = 0.01115579665134*q**3.11*(3.0*q**3 - 7.0*q**2 + 6.0)  &
+ 0.003923121822388*q**3.11*(3.0*q**3 - 6.0*q**2  &
+ 4.0) + 0.000881184569616*q**3.11*(33.0*q**3 - 84.0*q**2  &
+ 104.0) + 0.04009180980244*q**3.11 - 0.02004590490122*q**5.11  &
+ 0.007517214337957*q**6.11  &
- 0.1630499003148
    elseif (q <= 2.0) then
      df  = -0.003923121822388*q**3.11*(q - 2.0)**3 - 0.0004176230187754*q**3.11*(q  &
- 2.0)**2*(26.7126*q - 80.1378)  &
- 0.000881184569616*q**3.11*(11.0*q**3 - 84.0*q**2  &
+ 228.0*q - 208.0) + 0.08018361960487*q**3.11  &
- 0.06013771470365*q**4.11 + 0.02004590490122*q**5.11  &
- 0.002505738112652*q**6.11 - 0.1730728527654
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
      ddf  = 0.1666666666667*q**2.11*(3.0*q**3 - 6.0*q**2 + 4.0)
    elseif (q <= 2.0) then
      ddf  = -0.1666666666667*q**2.11*(q - 2.0)**3
    else
      ddf  = 0
    end if
  end subroutine
end module
