module cubic
  use const

  implicit none

  public :: kf, kdf, kddf, knorm, krad, kernelname

  private

    real :: knorm(3) = (/ 2./3., 10./(7. * pi), 1./(pi) /)
    real :: krad = 2.
    character (len=10) :: kernelname='cubic'

 contains

  subroutine kf(r, h, f)
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
      print *, 'something went wrong, q =', q
      stop
    end if
  end subroutine kf

  subroutine kdf(r, h, df)
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
      print *, 'something went wrong, q =', q
      stop
    end if
  end subroutine kdf

  subroutine kddf(r, h, ddf)
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
      print *, 'something went wrong, q =', q
      stop
    end if
  end subroutine kddf
end module cubic
