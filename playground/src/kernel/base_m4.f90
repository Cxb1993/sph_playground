module kernel_base
  use const

  implicit none

  public :: kf, kdf, kddf, knorm, krad, kernelname, initkernelbase, wCv, fwc, returnneibnum, cnarr, d2curnumb

  private

    real :: knorm(3) = [ 2./3., 10./(7. * pi), 1./(pi) ]
    real :: fwcl(3)  = [ 6., 6.3226, 6.66400]
    real :: cnarr(ecn_total)
    integer :: maxneibnum(3) = [35, 250, 150]

    real :: wCv, fwc, d2curnumb = -1.
    real :: krad = 2.
    integer :: dim, returnneibnum
    character (len=10) :: kernelname='M4'

 contains
   subroutine initkernelbase(d)
     integer, intent(in) :: d
     dim = d
     wCv = knorm(dim)
     fwc = fwcl(dim)
     returnneibnum = maxneibnum(dim)
     cnarr(ecn_hydro) = 0.1
     cnarr(ecn_d22nw) = 0.96
     cnarr(ecn_d2fab) = 0.15
     cnarr(ecn_d2n2w) = 0.2
   end subroutine initkernelbase

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
        error stop 'q is nan in kf'
      else
        error stop 'q is negative in kf'
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
        error stop 'q is nan in kdf'
      else
        error stop 'q is negative in kdf'
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
        error stop 'q is nan in kddf'
      else
        error stop 'q is negative in kddf'
      end if
    end if
  end subroutine kddf
end module kernel_base
