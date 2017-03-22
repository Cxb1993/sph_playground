module n2ext
  use const
  implicit none

  public :: n2f, n2df, n2ddf, n2C, n2R, n2kernelname!, calc_params

  private
    real, parameter :: n2R = 2.
    real :: n2C(3,2) = reshape(                       &
                       (/ 3./(sqrt(pi))/1.6,          &
                          30./(7.*pi) + 0.01,         &
                          50./(7.*pi**(3./2)) - 0.02, &
                          30./(7.*sqrt(pi)) - 0.0925, &
                          75./(14.*pi) + 0.052,       &
                          50./(7.*pi**(3./2)) + 0.1   &
                          /), (/3,2/))
    character (len=10) :: n2kernelname='mrunge'
 !    real, save :: f_err, df_err, ddf_err
 !
 contains
  subroutine n2f(r, h, f)
    real, intent(in)  :: r, h
    real, intent(out) :: f
    real              :: q

    q = r / h
    if (q >= n2R) then
      f  = 0.
    else if (q >= 0.) then
      f =  0.959609190185938 - 1.2000000000000006*q + 0.6004885122675785*q**3 -&
      0.30036638420068384*q**4 + 0.045073276840136775*q**5
      ! f = 0.9762537291293276 - 1.2147296804808763*q + 0.6165070397905309*q**3 -&
      ! 0.3215863301434455*q**4 + 0.05576150123907225*q**5 - 0.0018412100601094715*q**6
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
      df =  -1.2000000000000006 + 1.8014655368027355*q**2 - 1.2014655368027354*q**3 +&
      0.22536638420068386*q**4
      ! df = -1.2147296804808763 + 1.8495211193715928*q**2 - 1.286345320573782*q**3 + &
      ! 0.27880750619536127*q**4 - 0.011047260360656829*q**5
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
      ddf = 3.602931073605471*q - 3.604396610408206*q**2 + 0.9014655368027354*q**3
      ! ddf = 3.6990422387431856*q - 3.8590359617213466*q**2 + 1.115230024781445*q**3 - &
      ! 0.05523630180328414*q**4
    else
      print *, 'something went wrong, q =', q
      stop
    end if
  end subroutine n2ddf
end module n2ext
