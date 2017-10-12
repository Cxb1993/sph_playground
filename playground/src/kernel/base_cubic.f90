module base_kernel
  use const

  implicit none

  public :: kf, kdf, kddf, knorm, krad, kernelname, setdimbase, wCv, fwc

  private

    real :: knorm(3) = [ 2./3., 10./(7. * pi), 1./(pi) ]
    real :: fwcl(3)  = [ 6., 6.32261, 6.6666666 ]
    real :: krad = 2., wCv, fwc
    integer :: dim
    character (len=10) :: kernelname='cubic'

 contains
   subroutine setdimbase(d)
     integer, intent(in) :: d
     dim = d
     wCv = knorm(dim)
     fwc = fwcl(dim)
   end subroutine

  pure subroutine kf(q, f)
    real, intent(in)  :: q
    real, intent(out) :: f

    if (q >= 2.) then
      f  = 0.
    else if (q >= 1.) then
      f  = 0.25 * (2. - q)**3
    else if (q >= 0.) then
      f  = 0.25 * (2. - q)**3 - (1. - q)**3
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine kf

  pure subroutine kdf(q, df)
    real, intent(in)  :: q
    real, intent(out) :: df

    if (q >= 2.) then
      df = 0.
    else if (q >= 1.) then
      df = -0.75 * ((2. - q) ** 2)
    else if (q >= 0.) then
      df = -0.75 * (2. - q)**2 + 3 * (1. - q)**2
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine kdf

  pure subroutine kddf(q, ddf)
    real, intent(in)  :: q
    real, intent(out) :: ddf

    if (q >= 2.) then
      ddf = 0.
    else if (q >= 1.) then
      ddf = 1.5 * (2. - q)
    else if (q >= 0.) then
      ddf = 1.5 * (2. - q) - 6 * (1. - q)
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine kddf
end module
