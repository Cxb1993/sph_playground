module base_kernel
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
    elseif (q < 0.66666666) then
      ! f  = 15.0*(-1.5*q + 1.0)**5 - 6.0*(-1.5*q + 2.0)**5 + (-1.5*q  &
      !       + 3.0)**5
      f = 66. - 135.*q*q + 151.875*q*q*q*q - 75.9375*q*q*q*q*q
    elseif (q < 1.33333333) then
      ! f  = -6.0*(-1.5*q + 2.0)**5 + (-1.5*q + 3.0)**5
      f = 51. + 112.5*q - 472.5*q*q + 506.25*q*q*q - 227.813*q*q*q*q + 37.9688*q*q*q*q*q
    elseif (q < 2.0) then
      ! f  = (-1.5*q + 3.0)**5
      f = -7.59375*(-2. + q)**5
    else
      f  = 0.
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
    elseif (q < 0.66666666) then
      ! df  = q*(-379.6875*q*q*q + 607.5*q*q - 270.0)
      df = -270.*q + 607.5*q*q*q - 379.688*q*q*q*q
    elseif (q < 1.33333333) then
      ! df = -7.5*(1.5*q - 3.0)**4 + 45.0*(1.5*q - 2.0)**4
      df = 112.5 - 945.*q + 1518.75*q*q - 911.25*q*q*q + 189.844*q*q*q*q
    elseif (q < 2.0) then
      ! df  = -7.5*(1.5*q - 3.0)**4
      df  = -37.9687*(2. - q)**4
    else
      df  = 0.
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
    elseif (q < 0.66666666) then
      ! ddf  = -1518.75*q**3 + 1822.5*q**2 - 270.0
      ddf = -270. + 1822.5*q*q - 1518.75*q*q*q
    else if (q < 1.33333333) then
      ! ddf  = 759.375*q**3 - 2733.75*q**2 + 3037.5*q - 945.0
      ddf = -945. + 3037.5*q - 2733.75*q*q + 759.375*q*q*q
    else if (q < 2.0) then
      ! ddf  = -45.0*(1.5*q - 3.0)**3
      ddf = -151.875*(-2. + q)**3
    else
      ddf  = 0.
    end if
  end subroutine
end module
