module base_kernel
  use const
  implicit none

  public :: kf, kdf, kddf, knorm, krad, kernelname, setdimbase, wCv, fwc, returnneibnum

  private

    character (len=10) :: kernelname=' M6/3 '

    real :: knorm(3) = [ 1./120., 7./(478. * pi), 1./(120. * pi) ]
    real :: fwcl(3) = [4., 9.31505/2.25, 9.6429/2.25]
    integer :: maxneibnum(3) = [25, 570, 700]
    real :: krad = 3.

    real :: wCv, fwc
    integer :: dim, returnneibnum

 contains

  subroutine setdimbase(d)
    integer, intent(in) :: d

    dim = d
    wCv = knorm(dim)
    fwc = fwcl(dim)
    returnneibnum = maxneibnum(dim)
  end subroutine

  pure subroutine kf(q, f)
    real, intent(in)  :: q
    real, intent(out) :: f

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
  end subroutine

  pure subroutine kdf(q, df)
    real, intent(in)  :: q
    real, intent(out) :: df

    if (q >= 3.) then
      df = 0.
    else if (q >= 2.) then
      df = -5. * (3. - q)**4
    else if (q >= 1.) then
      df = -5. * (3. - q)**4 + 30. * (2. - q)**4
    else if (q >= 0.) then
      df = -5. * (3. - q)**4 + 30. * (2. - q)**4 - 75. * (1. - q)**4
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine

  pure subroutine kddf(q, ddf)
    real, intent(in)  :: q
    real, intent(out) :: ddf

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
  end subroutine
end module
