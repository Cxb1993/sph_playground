module n2ext
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase

  private

    real :: n2C(3) = (/ 6.66666682060241, 5.68446349992533, 4.71847639270427 /)
    character (len=10) :: kernelname = ' q*M4 '
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
      f  = -0.1666666666667*q + 0.001587301587302*q**4.0*(3.0*q**3  &
- 8.0*q**2 + 8.0) + 0.0005291005291005*q**4.0*(3.0*q**3  &
- 6.0*q**2 + 4.0) + 0.04074074074074*q**4.0  &
- 0.01746031746032*q**6.0 + 0.005555555555556*q**7.0  &
+ 0.147619047619
    elseif (q <= 2.0) then
      f  = -0.1269841269841*q - 0.0005291005291005*q**4.0*(q -  &
2.0)**3 - 6.613756613757e-5*q**4.0*(q - 2.0)**2*(24.0*q  &
- 96.0) + 0.03174603174603*q**4.0 - 0.01587301587302*q**5.0  &
+ 0.003968253968254*q**6.0 - 0.0003968253968254*q**7.0  &
- 6.613756613757e-5*(2.0*q -  &
4.0)*(188.0*q**4.0 - 80.0*q**5.0 + 11.0*q**6.0 +  &
384.0) + 0.05079365079365
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
      df  = 0.04444444444444*q**3.0 + 0.01111111111111*q**3.0*(3.0*q**3  &
- 7.0*q**2 + 6.0) + 0.003703703703704*q**3.0*(3.0*q**3  &
- 6.0*q**2 + 4.0) + 0.0009259259259259*q**3.0*(33.0*q**3  &
- 84.0*q**2 + 104.0) - 0.02222222222222*q**5.0  &
+ 0.008333333333333*q**6.0 - 0.1666666666667
    elseif (q <= 2.0) then
      df  = 0.08888888888889*q**3.0 - 0.003703703703704*q**3.0*(q  &
- 2.0)**3 - 0.000462962962963*q**3.0*(q - 2.0)**2*(24.0*q  &
- 72.0) - 0.0009259259259259*q**3.0*(11.0*q**3  &
- 84.0*q**2 + 228.0*q - 208.0) - 0.06666666666667*q**4.0  &
+ 0.02222222222222*q**5.0 - 0.002777777777778*q**6.0  &
- 0.1777777777778
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
      ddf  = 0.1666666666667*q**2.0*(3.0*q**3 - 6.0*q**2 + 4.0)
    elseif (q <= 2.0) then
      ddf  = -0.1666666666667*q**2.0*(q - 2.0)**3
    else
      ddf  = 0
    end if
  end subroutine
end module n2ext
