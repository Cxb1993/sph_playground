module kernel_base
  use const
  implicit none

  public :: kf, kdf, kddf, knorm, krad, kernelname, setdimbase, wCv, fwc, returnneibnum

  private

    character (len=10) :: kernelname='M8'

    real :: knorm(3) = [1., 2268./(1487.*pi), 3./(4.*pi)]
    real :: fwcl(3) = [3., 3.07775455, 3.15791438]
    integer :: maxneibnum(3) = [100, 200, 400]
    real :: krad = 4.

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

    if (q >= 4.) then
      f = 0.
    else if (q >= 3.) then
      f = -(-4 + q)**7/5040.
    else if (q >= 2.) then
      f = (-1112 + 12152*q - 19320*q**2 + 13720*q**3 - 5320*q**4 + 1176*q**5 - 140*q**6 + 7*q**7)/5040.
    else if (q >= 1.) then
      f = (2472 - 392*q - 504*q**2 - 1960*q**3 + 2520*q**4 - 1176*q**5 + 252*q**6 - 21*q**7)/5040.
    else if (q >= 0.) then
      f = 0.4793650793650794 - q**2/3. + q**4/9. - q**6/36. + q**7/144.
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

    if (q >= 4.) then
      df = 0.
    else if (q >= 3.) then
      df = -(-4 + q)**6/720.
    else if (q >= 2.) then
      df = (12152 - 38640*q + 41160*q**2 - 21280*q**3 + 5880*q**4 - 840*q**5 + 49*q**6)/5040.
    else if (q >= 1.) then
      df = (-392 - 1008*q - 5880*q**2 + 10080*q**3 - 5880*q**4 + 1512*q**5 - 147*q**6)/5040.
    else if (q >= 0.) then
      df = (-2*q)/3. + (4*q**3)/9. - q**5/6. + (7*q**6)/144.
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

    if (q >= 4.) then
      ddf = 0.
    else if (q >= 3.) then
      ddf = -(-4 + q)**5/120.
    else if (q >= 2.) then
      ddf = (-38640 + 82320*q - 63840*q**2 + 23520*q**3 - 4200*q**4 + 294*q**5)/5040.
    else if (q >= 1.) then
      ddf = (-1008 - 11760*q + 30240*q**2 - 23520*q**3 + 7560*q**4 - 882*q**5)/5040.
    else if (q >= 0.) then
      ddf = -0.6666666666666666 + (4*q**2)/3. - (5*q**4)/6. + (7*q**5)/24.
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine
end module kernel_base
