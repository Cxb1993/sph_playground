module n2movedgaus
  use const
  implicit none

  public :: n2f, n2df, n2ddf, n2C, n2R, n2kernelname!, calc_params

  private
    real, parameter :: n2R = 2.
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
    ! real :: knorm(3) = (/ 1./(2*sqrt(pi) - 0.17), &
    !                       1./(pi - 0.1214),     &
    !                       1./(pi**(3./2) - 0.375) /)
    ! R = 2 Own
    ! real :: knorm(3) = (/ 1./(sqrt(pi) - 0.835), &
    !                       1./(pi - 0.115),     &
    !                       1./(pi**(3./2) - .29) /)

    real :: n2C(3,2) = reshape(                       &
                       (/ 3./(sqrt(pi)) - 0.075,      &
                          30./(7.*pi) + 0.01,         &
                          50./(7.*pi**(3./2)) - 0.024,&
                          30./(7.*sqrt(pi)) - 0.0925, &
                          75./(14.*pi) + 0.052,       &
                          50./(7.*pi**(3./2)) + 0.1   &
                          /), (/3,2/))
    character (len=10) :: n2kernelname='mgauss'
 !    real, save :: f_err, df_err, ddf_err
 !
 contains
 !  subroutine calc_params
 !    real fa, fr, q
 !    ! Where we want to truncate
 !    q = krad
 !    ! What we have at the truncation point
 !    fr = exp(-q**2) * (4 * q**2 - 2)
 !    ! What we have to have at the real truncation point ie 'inf'
 !    fa = 0
 !    ! Make it the same
 !    ddf_err = fr - fa
 !    ! But with this we change the undercurv volume, so we need to change it back
 !    ! Page 26 + something strange
 !    fr = (-2*exp(-q**2) - ddf_err/2/krad*q)*q
 !    fa = 0
 !    df_err  = fr - fa
 !
 !    fr = exp(-q**2) - (q**2*(ddf_err - 3*df_err))/(6*krad)
 !    fa = 0
 !    f_err = fr - fa
 !
 !    ! print *, ddf_err, df_err, f_err
 !  end subroutine calc_params
  subroutine n2f(r, h, f)
    real, intent(in)  :: r, h
    real, intent(out) :: f
    real              :: q

    q = r / h
    if (q >= n2R) then
      f  = 0.
    else if (q >= 0.) then
      ! f  = exp(-q**2)
      ! e^(-q^2) - (q^2 (e_1 q - 3 e_2))/(6 R)
      ! f = exp(-q**2) - (q**2*(ddf_err - 3*df_err))/(6*krad) - q/krad*f_err
      ! f = 0.259118 + 0.0131749*exp(6.*q - 4.*q**2) - 0.335603*q - 0.0526996*q**2 + 0.0086224*q**3 +&
      ! 0.332335*Erf(1.5 - 2.*q) - 0.443113*q*Erf(1.5 - 2.*q)

      f = 0.259118 + 0.00966309*2.71828**((6.4 - 4.*q)*q) + &
          q*(-0.362353 + (-0.0386524 + 0.00617947*q)*q) + &
          (0.354491 - 0.443113*q)*erf(1.6 - 2.*q)
    else
      print *, 'something went wrong, q =', q
      stop
    end if
  end subroutine n2f

  subroutine n2df(r, h, df)
    real, intent(in)  :: r, h
    real, intent(out) :: df
    real              :: q

    q = r / h
    if (q >= n2R) then
      df = 0.
    else if (q >= 0.) then
      ! all this also multiplied on 'q'
      ! df = -2 * exp(-q**2)
      ! (-2 e^(-q^2) - ddf_err/2/R*q)*q - df_err*q/3
      ! -2 e^(-q^2) - ddf_err/2/R*q - df_err/3  | :q
      ! df = -2 *exp(-q**2) - ddf_err/2/krad*q - df_err/krad
      ! Gaus bump on 1 linear
      ! df = 2./exp(4.) - q/exp(4.) - 1./4.*sqrt(pi)*erf(2.) + 1./4.*sqrt(pi) * erf(2.*(q - 1.))
      ! Gaus bump on 1 quadric
      ! df = -0.433714 - 0.000228945*(-2. + q)**5 - 0.000228945*q**5 - 0.443113*erf(2. - 2.*q)
      ! WOOOW
      df = -0.362353 + (-0.0773047 + 0.0185384*q)*q - 0.443113*erf(1.6 - 2.*q)
    else
      print *, 'something went wrong, q =', q
      stop
    end if
  end subroutine n2df

  subroutine n2ddf(r, h, ddf)
    real, intent(in)  :: r, h
    real, intent(out) :: ddf
    real              :: q!, fr, fa

    q = r / h
    if (q >= n2R) then
      ddf = 0.
    else if (q >= 0.) then
      ! ddf = exp(-q**2) * (4 * q**2 - 2)
      ! ddf = exp(-((q-krad/2)*krad)**2) - 1/exp(4.)
      ! Gaus bump is moved to 1 with linear correction page 27
      ! ddf = exp(-((q - 1.)*2)**2) - 1./exp(4.) * q/krad - 1./exp(4.) * (1 - q/krad)
      ! Gaus bump is moved to 1 with quasric corrections page 27
      ! ddf = exp(-4.*(-1. + q)**2) - (1. - q/2.)**4/exp(4.) - q**4/(16*exp(4.))
      ! WOOOW
      ddf = exp(-4*(-0.8 + q)**2) - 0.0773047*(1. - q/2.) - 0.00157556*q
    else
      print *, 'something went wrong, q =', q
      stop
    end if
  end subroutine n2ddf
end module n2movedgaus
