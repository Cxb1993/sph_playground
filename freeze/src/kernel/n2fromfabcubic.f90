module n2fromfabcubic
  use const
  implicit none

  public :: n2f, n2df, n2ddf, n2R, n2Name, setdimbase, n2Cv


  private
    real, parameter :: n2C(3) = (/  2./3., 10./(7. * pi), 1./(pi) /)

    character (len=10) :: n2Name = ' genesis '
    real               :: n2R = 2.0, n2Cv
    integer            :: dim

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
    if (dim == 1) then
      ! (0, q > 2),
      ! (0.25*q**3 - 3.0*q**2 + 6.0*q*log(q) - 1.15888308335967*q + 4.0, q > 1),
      ! (-0.75*q**3 + 3.0*q**2 - 4.15888308335967*q + 2.0, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2f'
        end if
      elseif ( q < 1.0 ) then
        f = -0.75*q**3 + 3.0*q**2 - 4.15888308335967*q + 2.0
      elseif ( q < 2.0 ) then
        f = 0.25*q**3 - 3.0*q**2 + 6.0*q*log(q) - 1.15888308335967*q + 4.0
      else
        f = .0
      end if
    elseif ( dim == 2 ) then
      ! (0, q > 2)
      ! (0.1666666*q**3 - 1.5*q**2 + 6.0*q - 3.9999984*log(q) - 4.56074518679571, q > 1)
      ! (-0.5*q**3 + 1.5*q**2 - 1.9999986*log(q) - 0.894078586795709, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2f'
        end if
      elseif ( q < 1.0 ) then
        f = -0.5*q**3 + 1.5*q**2 - 2.0*log(q) - 0.894078586795709
      elseif ( q < 2.0 ) then
        f  = 1./6.*q**3 - 1.5*q**2 + 6.0*q - 4.*log(q) - 4.56074518679571
      else
        f = .0
      end if
    elseif ( dim == 3 ) then
      ! (0, q > 2)
      ! (0.125*q**3 - q**2 + 3*q - 4.0 + 2.0/q, q > 1)
      ! (-0.375*q**3 + q**2 - 2.0 + 1.5/q, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2f'
        end if
      elseif ( q < 1.0 ) then
        f = -0.375*q**3 + q**2 - 2.0 + 1.5/q
      elseif ( q < 2.0 ) then
        f  = 0.125*q**3 - q**2 + 3*q - 4.0 + 2.0/q
      else
        f = .0
      end if
    end if
  end subroutine

  ! pure subroutine n2df(r, h, df)
  pure subroutine n2df(r, h, df)
    real, intent(in)  :: r, h
    real, intent(out) :: df
    real              :: q

    q = r / h
    if (dim == 1) then
      ! (0, q > 2)
      ! (0.75*q**2 - 6.0*q + 6.0*log(q) + 4.84111691664033, q > 1)
      ! (-2.25*q**2 + 6.0*q - 4.15888308335967, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2df'
        end if
      elseif ( q < 1.0 ) then
        df = -2.25*q**2 + 6.0*q - 4.15888308335967
      elseif ( q < 2.0 ) then
        df = 0.75*q**2 - 6.0*q + 6.0*log(q) + 4.84111691664033
      else
        df = .0
      end if
    elseif ( dim == 2 ) then
      ! (0, q > 2)
      ! (0.4999998*q**2 - 3.0*q + 6.0 - 3.9999984/q, q > 1)
      ! (-1.5*q**2 + 3.0*q - 1.9999986/q, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2df'
        end if
      elseif ( q < 1.0 ) then
        df = -1.5*q**2 + 3.0*q - 2.0/q
      elseif ( q < 2.0 ) then
        df = 0.5*q**2 - 3.0*q + 6.0 - 4.0/q
      else
        df = .0
      end if
    elseif ( dim == 3 ) then
      ! (0, q > 2)
      ! (0.375*q**2 - 2*q + 3 - 2.0/q**2, q > 1)
      ! (-1.125*q**2 + 2*q - 1.5/q**2, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2df'
        end if
      elseif ( q < 1.0 ) then
        df = -1.125*q**2 + 2*q - 1.5/q**2
      elseif ( q < 2.0 ) then
        df = 0.375*q**2 - 2*q + 3 - 2.0/q**2
      else
        df = .0
      end if
    end if
    df = df / q
  end subroutine

  ! pure subroutine n2ddf(r, h, ddf)
  pure subroutine n2ddf(r, h, ddf)
    real, intent(in)  :: r, h
    real, intent(out) :: ddf
    real              :: q

    q = r / h
    if (dim == 1) then
      ! (0, q > 2)
      ! (1.5*q - 6.0 + 6.0/q, q > 1)
      ! (-4.5*q + 6.0, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2ddf'
        end if
      elseif ( q < 1.0 ) then
        ddf = -4.5*q + 6.0
      elseif ( q < 2.0 ) then
        ddf = 1.5*q - 6.0 + 6.0/q
      else
        ddf = .0
      end if
    elseif ( dim == 2 ) then
      !  (0, q > 2),
      !  (1.0*q - 3.0 + 4.0/q**2, q > 1),
      !  (-3.0*q + 3.0 + 2.0/q**2, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2ddf'
        end if
      elseif ( q < 1.0 ) then
        ddf = -3.0*q + 3.0 + 2.0/q**2
      elseif ( q < 2.0 ) then
        ddf = 1.0*q - 3.0 + 4.0/q**2
      else
        ddf = .0
      end if
    elseif ( dim == 3 ) then
      ! (0, q > 2)
      ! (0.75*q - 2 + 4.0/q**3, q > 1)
      ! (-2.25*q + 2 + 3.0/q**3, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2ddf'
        end if
      elseif ( q < 1.0 ) then
        ddf = -2.25*q + 2 + 3.0/q**3
      elseif ( q < 2.0 ) then
        ddf = 0.75*q - 2 + 4.0/q**3
      else
        ddf = .0
      end if
    end if
  end subroutine
end module
