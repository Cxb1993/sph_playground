module gaus
  implicit none

  public :: kf, kdf, kddf, knorm, krad, calc_params

  private
    real, parameter :: pi = 4.*atan(1.)
    real, parameter :: krad = 2.
    ! double Q
    ! real :: knorm(3) = (/ 8./sqrt(pi),  &
    !                       16./pi,        &
    !                       32./pi**(3./2) &
    !                     /)
    ! real :: knorm(3) = (/ 1./sqrt(pi),  &
    !                       1./pi,        &
    !                       1./pi**(3./2) &
    !                     /)
    ! R = 3
    ! real :: knorm(3) = (/ 1./(sqrt(pi) - 0.045), &
    !                       1./(pi - 0.115),     &
    !                       1./(pi**(3./2) - .29) /)
    ! R = 2
    real :: knorm(3) = (/ 1./(sqrt(pi) - 1.), &
                          1./(pi - 0.115),     &
                          1./(pi**(3./2) - .29) /)

    real, save :: f_err, df_err, ddf_err

 contains
  subroutine calc_params
    real fa, fr, q
    ! Where we want to truncate
    q = krad
    ! What we have at the truncation point
    fr = exp(-q**2) * (4 * q**2 - 2)
    ! What we have to have at the real truncation point ie 'inf'
    fa = 0
    ! Make it the same
    ddf_err = fr - fa
    ! But with this we change the undercurv volume, so we need to change it back
    ! Page 26 + something strange
    fr = -2 * q * exp(-q**2) - ddf_err * q
    fa = 0
    df_err  = fr - fa

    fr = exp(-q**2) - ddf_err * q - df_err * log(q)
    fa = 0
    f_err = fr - fa

    ! print *, ddf_err, df_err, f_err
  end subroutine calc_params

  subroutine kf(r, h, f)
    real, intent(in)  :: r, h
    real, intent(out) :: f
    real              :: q

    q = r / h
    if (q >= krad) then
      f  = 0.
    else if (q >= 0.) then
      ! f  = exp(-q**2)
      f = exp(-q**2) - ddf_err * q - df_err * log(q) - f_err
    else
      print *, 'something went wrong, q =', q
      stop
    end if
  end subroutine kf

  subroutine kdf(r, h, df)
    real, intent(in)  :: r, h
    real, intent(out) :: df
    real              :: q

    q = r / h
    if (q >= krad) then
      df = 0.
    else if (q >= 0.) then
      ! all this also multiplied on 'q'
      ! df = -2 * exp(-q**2)
      df = -2 *exp(-q**2) - ddf_err - df_err/q
    else
      print *, 'something went wrong, q =', q
      stop
    end if
  end subroutine kdf

  subroutine kddf(r, h, ddf)
    real, intent(in)  :: r, h
    real, intent(out) :: ddf
    real              :: q, fr, fa

    q = r / h
    if (q >= krad) then
      ddf = 0.
    else if (q >= 0.) then
      ! ddf = exp(-q**2) * (4 * q**2 - 2)
      ddf = exp(-q**2) * (4 * q**2 - 2) - ddf_err
    else
      print *, 'something went wrong, q =', q
      stop
    end if
  end subroutine kddf
end module gaus
