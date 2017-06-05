module n2fromWcubic
  use const
  implicit none

  public :: kf, kdf, kddf, krad, kernelname, setdimbase, wCv


  private
    real, parameter :: n2C(3) = (/ 4.00000006013221, 2.87495346287778, 2.12206600827867 /)

    character (len=10) :: kernelname = ' genesis '
    real               :: krad = 2.0, wCv
    integer            :: dim

  contains

  subroutine setdimbase(d)
    integer, intent(in) :: d

    dim = d
    wCv = n2C(dim)
  end subroutine

  pure subroutine kf(q, f)
    real, intent(in)  :: q
    real, intent(out) :: f

    if (dim == 1) then
      ! (0, q > 2)
      ! (-0.0125*q**5 + 0.125*q**4 - 0.5*q**3 + 1.0*q**2 - 1.0*q + 0.4, q > 1)
      ! (0.0375*q**5 - 0.125*q**4 + 0.5*q**2 - 0.75*q + 0.35, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2f'
        end if
      elseif ( q < 1.0 ) then
        f = 0.0375*q**5 - 0.125*q**4 + 0.5*q**2 - 0.75*q + 0.35
      elseif ( q < 2.0 ) then
        f = -0.0125*q**5 + 0.125*q**4 - 0.5*q**3 + 1.0*q**2 - 1.0*q + 0.4
      else
        f = .0
      end if
    elseif ( dim == 2 ) then
      ! (0, q > 2)
      ! (-0.01*q**5 + 0.09375*q**4 - 0.333333333333333*q**3 + 0.5*q**2 - 0.4*log(q) - 0.236074461109355, q > 1)
      ! (0.03*q**5 - 0.09375*q**4 + 0.25*q**2 - 0.35*log(q) - 0.171907794442689, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2f'
        end if
      elseif ( q < 1.0 ) then
        f = 0.03*q**5 - 0.09375*q**4 + 0.25*q**2 - 0.35*log(q) - 0.171907794442689
      elseif ( q < 2.0 ) then
        f = -0.01*q**5 + 0.09375*q**4 - 0.333333333333333*q**3 + 0.5*q**2 - 0.4*log(q) - 0.236074461109355
      else
        f = .0
      end if
    elseif ( dim == 3 ) then
      ! (0, q > 2)
      ! (-0.00833333333333333*q**5 + 0.075*q**4 - 0.25*q**3 + 0.333333333333333*q**2 - 0.4 + 0.266666666666667/q, q > 1)
      ! (0.025*q**5 - 0.075*q**4 + 0.166666666666667*q**2 - 0.35 + 0.25/q, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2f'
        end if
      elseif ( q < 1.0 ) then
        f = 0.025*q**5 - 0.075*q**4 + 0.166666666666667*q**2 - 0.35 + 0.25/q
      elseif ( q < 2.0 ) then
        f = -0.00833333333333333*q**5 + 0.075*q**4 - 0.25*q**3 + 0.333333333333333*q**2 - 0.4 + 0.266666666666667/q
      else
        f = .0
      end if
    end if
  end subroutine

  pure subroutine kdf(q, df)
    real, intent(in)  :: q
    real, intent(out) :: df

    if (dim == 1) then
      ! (0, q > 2)
      ! (-0.0625*q**4 + 0.5*q**3 - 1.5*q**2 + 2.0*q - 1.0, q > 1)
      ! (0.1875*q**4 - 0.5*q**3 + 1.0*q - 0.75, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2df'
        end if
      elseif ( q < 1.0 ) then
        df = 0.1875*q**4 - 0.5*q**3 + 1.0*q - 0.75
      elseif ( q < 2.0 ) then
        df = -0.0625*q**4 + 0.5*q**3 - 1.5*q**2 + 2.0*q - 1.0
      else
        df = .0
      end if
    elseif ( dim == 2 ) then
      ! (0, q > 2)
      ! (-0.05*q**4 + 0.375*q**3 - 1.0*q**2 + 1.0*q - 0.4/q, q > 1)
      ! (0.15*q**4 - 0.375*q**3 + 0.5*q - 0.35/q, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2df'
        end if
      elseif ( q < 1.0 ) then
        df = 0.15*q**4 - 0.375*q**3 + 0.5*q - 0.35/q
      elseif ( q < 2.0 ) then
        df = -0.05*q**4 + 0.375*q**3 - 1.0*q**2 + 1.0*q - 0.4/q
      else
        df = .0
      end if
    elseif ( dim == 3 ) then
      ! (0, q > 2)
      ! (-0.0416666666666667*q**4 + 0.3*q**3 - 0.75*q**2 + 0.666666666666667*q - 0.266666666666667/q**2, q > 1)
      ! (0.125*q**4 - 0.3*q**3 + 0.333333333333333*q - 0.25/q**2, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2df'
        end if
      elseif ( q < 1.0 ) then
        df = 0.125*q**4 - 0.3*q**3 + 0.333333333333333*q - 0.25/q**2
      elseif ( q < 2.0 ) then
        df = -0.0416666666666667*q**4 + 0.3*q**3 - 0.75*q**2 + 0.666666666666667*q - 0.266666666666667/q**2
      else
        df = .0
      end if
    end if
  end subroutine

  pure subroutine kddf(q, ddf)
    real, intent(in)  :: q
    real, intent(out) :: ddf

    if (dim == 1) then
      ! (0, q > 2)
      ! (-0.25*q**3 + 1.5*q**2 - 3.0*q + 2.0, q > 1)
      ! (0.75*q**3 - 1.5*q**2 + 1.0, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2ddf'
        end if
      elseif ( q < 1.0 ) then
        ddf = 0.75*q**3 - 1.5*q**2 + 1.0
      elseif ( q < 2.0 ) then
        ddf = -0.25*q**3 + 1.5*q**2 - 3.0*q + 2.0
      else
        ddf = .0
      end if
    elseif ( dim == 2 ) then
      ! (0, q > 2)
      ! (-0.2*q**3 + 1.125*q**2 - 2.0*q + 1.0 + 0.4/q**2, q > 1)
      ! (0.6*q**3 - 1.125*q**2 + 0.5 + 0.35/q**2, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2ddf'
        end if
      elseif ( q < 1.0 ) then
        ddf = 0.6*q**3 - 1.125*q**2 + 0.5 + 0.35/q**2
      elseif ( q < 2.0 ) then
        ddf = -0.2*q**3 + 1.125*q**2 - 2.0*q + 1.0 + 0.4/q**2
      else
        ddf = .0
      end if
    elseif ( dim == 3 ) then
      ! (0, q > 2)
      ! (-0.166666666666667*q**3 + 0.9*q**2 - 1.5*q + 0.666666666666667 + 0.533333333333333/q**3, q > 1)
      ! (0.5*q**3 - 0.9*q**2 + 0.333333333333333 + 0.5/q**3, q > 0)
      if (q <= 0.0) then
        if (isnan(q)) then
          error stop 'q is nan'
        else
          error stop 'q is negative n2ddf'
        end if
      elseif ( q < 1.0 ) then
        ddf = 0.5*q**3 - 0.9*q**2 + 0.333333333333333 + 0.5/q**3
      elseif ( q < 2.0 ) then
        ddf = -0.166666666666667*q**3 + 0.9*q**2 - 1.5*q + 0.666666666666667 + 0.533333333333333/q**3
      else
        ddf = .0
      end if
    end if
  end subroutine
end module
