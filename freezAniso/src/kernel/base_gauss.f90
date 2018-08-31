module base_kernel
  use const

  implicit none

  public :: kf, kdf, kddf, knorm, krad, kernelname, setdimbase, wCv, fwc, returnneibnum

  private

    real :: knorm(3) = [ 1./sqrt(pi), 1./pi, 1./pi**(1.5) ]
    real :: fwcl(3) = [4., 4., 4.]
    integer :: maxneibnum(3) = [1, 1, 1]

    real :: krad = 2., wCv, fwc
    integer :: dim, returnneibnum
    character (len=10) :: kernelname='gauss'

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

    if (q >= 0.) then
      f  = exp(-q * q)
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

    if (q >= 0.) then
      df = -2. * exp(-q * q) * q
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

    if (q >= 0.) then
      ddf = -2 * exp(-q*q) + 4 * exp(-q*q) * q * q
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine kddf
end module
