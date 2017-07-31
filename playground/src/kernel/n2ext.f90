module n2ext
  use const
  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase, fwc

  private

    real :: n2C(3) = (/ 0.0124999997029078, 0.0104882435025047, 0.00895246473927735 /)
    real :: fwcl(3) = [9., 9.3151, 9.6429]
    character (len=10) :: kernelname = ' zipM6 '
    real :: krad = 2.0, wCv, fwc
    integer :: dim
  contains


    ! F1[q_] = M6q[q, 1];
    ! F2[q_] = -DDM6[q, 1];
    ! F3[q_] = M6q[q, 1];

  subroutine setdimbase(d)
    integer, intent(in) :: d
    dim = d
    wCv = n2C(dim)
    fwc = fwcl(dim)
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
    elseif (1.5*q < 1.0) then
      f  = 15.0*(-1.5*q + 1.0)**5 - 6.0*(-1.5*q + 2.0)**5 + (-1.5*q  &
+ 3.0)**5
    elseif (1.5*q < 2.0) then
      f  = -6.0*(-1.5*q + 2.0)**5 + (-1.5*q + 3.0)**5
    elseif (1.5*q < 3.0) then
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
    elseif (1.5*q < 1.0) then
      df  = q*(-379.6875*q**3 + 607.5*q**2 - 270.0)
    elseif (1.5*q < 2.0) then
      df  = -7.5*(1.5*q - 3.0)**4 + 45.0*(1.5*q - 2.0)**4
    elseif (1.5*q < 3.0) then
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
    elseif (1.5*q < 1.0) then
      ddf  = -1518.75*q**3 + 1822.5*q**2 - 270.0
    else if (1.5*q < 2.0) then
      ddf  = 759.375*q**3 - 2733.75*q**2 + 3037.5*q - 945.0
    else if (1.5*q < 3.0) then
      ddf  = -45.0*(1.5*q - 3.0)**3
    else
      ddf  = 0
    end if
  end subroutine
end module
