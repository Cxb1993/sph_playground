module n2ext
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase

  private

    real :: n2C(3) = [0.0124999997029078, 0.0104882435025047, 0.00895246473927735]
    character (len=10) :: kernelname = '  mmq0zipm6  '
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
    elseif (q < 0.666666666666667) then
      f  = -4.5*q**2 - 69.0*q + 74.0
    elseif (q < 1.33333333333333) then
      f  = -6.0*(-1.5*q + 2.0)**5 + (-1.5*q + 3.0)**5
    elseif (q < 2.0) then
      f  = (-1.5*q + 3.0)**5
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
    elseif (q < 0.666666666666667) then
      df  = -9.0*q - 69.0
    elseif (q < 1.33333333333333) then
      df  = -7.5*(1.5*q - 3.0)**4 + 45.0*(1.5*q - 2.0)**4
    elseif (q < 2.0) then
      df  = -7.5*(1.5*q - 3.0)**4
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
    elseif (q < 0.666666666666667) then
      ddf  = -9.000000000000
    else if (q < 1.33333333333333) then
      ddf  = 759.375*q**3 - 2733.75*q**2 + 3037.5*q - 945.0
    else if (q < 2.0) then
      ddf  = -45.0*(1.5*q - 3.0)**3
    else
      ddf  = 0
    end if
  end subroutine
end module
