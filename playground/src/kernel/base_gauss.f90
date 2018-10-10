module kernel_base
  use const

  implicit none

  public :: kf, kdf, kddf, knorm, krad, kernelname, initkernelbase, wCv, fwc, returnneibnum, cnarr, d2curnumb


  private

    real :: knorm(3) = [ 1./sqrt(pi), 1./pi, 1./pi**(1.5) ]
    real :: fwcl(3) = [4., 4., 4.]
    real :: cnarr(ecn_total)
    integer :: maxneibnum(3) = [100, 500, 1000]

    real :: wCv, fwc, d2curnumb = -1.
    real :: krad = 10.
    integer :: dim, returnneibnum
    character (len=10) :: kernelname='Gauss'

 contains

  subroutine initkernelbase(d)
    integer, intent(in) :: d

    dim = d
    wCv = knorm(dim)
    fwc = fwcl(dim)
    returnneibnum = maxneibnum(dim)
    cnarr(ecn_hydro) = 0.2
    cnarr(ecn_d22nw) = 0.6
    cnarr(ecn_d2fab) = 0.15
    cnarr(ecn_d2n2w) = 0.2
  end subroutine initkernelbase

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
end module kernel_base
