module n2ext
  use const
  implicit none

  public :: n2f, n2df, n2ddf, n2Cv, n2R, n2Name, setdimbase

  private

    real :: n2C(3) = (/ 4.00000009317521, 4.31258589868469, 4.24413248021448 /)
    character (len=10) :: n2Name=' genesis '
    real :: n2R = 2.0, n2Cv
    integer :: dim
  contains

  subroutine setdimbase(d)
    integer, intent(in) :: d
    dim = d
    n2Cv = n2C(dim)
  end subroutine

  pure subroutine n2f(r, h, f)
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
    elseif (q < 1.0) then
      f  = 0.0375*q**5 - 0.125*q**4 + 0.5*q**2 - 0.75*q + 0.35
    elseif (q < 2.0) then
      f  = -0.0125*q**5 + 0.125*q**4 - 0.5*q**3 + 1.0*q**2 - 1.0*q + 0.4
    else
      f  = 0
    end if
  end subroutine

  pure subroutine n2df(r, h, df)
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
    elseif (q < 1.0) then
      df  = (0.1875*q**4 - 0.5*q**3 + 1.0*q - 0.75)/q
    elseif (q < 2.0) then
      df  = (-0.0625*q**4 + 0.5*q**3 - 1.5*q**2 + 2.0*q - 1.0)/q
    else
      df  = 0
    end if
  end subroutine n2df

  pure subroutine n2ddf(r, h, ddf)
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
    elseif (q < 1.0) then
      ddf  = 0.75*q**3 - 1.5*q**2 + 1.0
    elseif (q < 2.0) then
      ddf  = -0.25*q**3 + 1.5*q**2 - 3.0*q + 2.0
    else
      ddf  = 0
    end if
  end subroutine n2ddf
end module n2ext
