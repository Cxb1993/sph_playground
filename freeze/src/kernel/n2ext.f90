module n2ext
  use const
  implicit none

  public :: n2f, n2df, n2ddf, n2Cv, n2R, n2Name, setdimbase


  private
    real :: n2C(3) = (/ 2./3., 10./(7. * pi), 1./(pi) /)

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
    if (q >= 2.) then
      f  = 0.
    else if (q >= 1.) then
      f  = 0.25 * (2. - q)**3
    else if (q >= 0.) then
      f  = 0.25 * (2. - q)**3 - (1. - q)**3
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine n2f

  pure subroutine n2df(r, h, df)
    real, intent(in)  :: r, h
    real, intent(out) :: df
    real              :: q

    q = r / h
    if (q >= 2.) then
      df = 0.
    else if (q >= 1.) then
      df = (- 0.75 * ((2. - q) ** 2)) / q
    else if (q >= 0.) then
      df = -3. + 2.25 * q
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine n2df

  pure subroutine n2ddf(r, h, ddf)
    real, intent(in)  :: r, h
    real, intent(out) :: ddf
    real              :: q

    q = r / h
    if (q >= 2.) then
      ddf = 0.
    else if (q >= 1.) then
      ddf = 3. - 1.5 * q
    else if (q >= 0.) then
      ddf = -3. + 4.5 * q
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine n2ddf
end module n2ext
