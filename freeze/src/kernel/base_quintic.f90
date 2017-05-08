module quintic
  use const
  implicit none

  public :: kf, kdf, kddf, knorm, krad, kernelname

  private

    real :: knorm(3) = (/ 1./120., 7./(478. * pi), 1./(120. * pi) /)
    real :: krad = 3.
    character (len=10) :: kernelname='quintic'

 contains

  pure subroutine kf(r, h, f)
    real, intent(in)  :: r, h
    real, intent(out) :: f
    real              :: q

    q = r / h
    if (q >= 3.) then
      f = 0.
    else if (q >= 2.) then
      f = (3. - q)**5
    else if (q >= 1.) then
      f = (3. - q)**5 - 6. * (2. - q)**5
    else if (q >= 0.) then
      f = (3. - q)**5 - 6. * (2. - q)**5 + 15. * (1. - q)**5
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine kf

  pure subroutine kdf(r, h, df)
    real, intent(in)  :: r, h
    real, intent(out) :: df
    real              :: q

    q = r / h
    if (q >= 3.) then
      df = 0.
    else if (q >= 2.) then
      df = -5. * (3. - q)**4 / q
    else if (q >= 1.) then
      df = (-5. * (3. - q)**4 + 30. * (2. - q)**4) / q
    else if (q > 0.) then
      df = (-5. * (3. - q)**4 + 30. * (2. - q)**4 - 75. * (1. - q)**4) / q
    else if (q == 0.) then
      df = 0.
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine kdf

  pure subroutine kddf(r, h, ddf)
    real, intent(in)  :: r, h
    real, intent(out) :: ddf
    real              :: q

    q = r / h
    if (q >= 3.) then
      ddf = 0.
    else if (q >= 2.) then
      ddf = 20. * (3. - q)**3
    else if (q >= 1.) then
      ddf = 20. * (3. - q)**3 - 120. * (2. - q)**3
    else if (q >= 0.) then
      ddf = 20. * (3. - q)**3 - 120. * (2. - q)**3 + 300. * (1. - q)**3
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine kddf
end module quintic
