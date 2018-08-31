module base_kernel
  use const

  implicit none

  public :: kf, kdf, kddf, knorm, krad, kernelname, setdimbase, wCv, fwc

  private

    real :: knorm(3) = [ 0., 1./0.0804323, 0. ]
    real :: fwcl(3)  = [ 0., 0., 0. ]
    real :: krad = 2., wCv, fwc
    integer :: dim
    character (len=10) :: kernelname='uspeh'

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
      f = 0.2365842017602397 + q**2*(-0.4693287500000001 +&
          q*(0.3446930555555555 + (-0.098861875 +&
          0.010000000000000002*q)*q)) + 0.20915666666666688*Log(q)
    else if (q > eps0) then
      f = 0.17241753509357294 + q**2*(-0.21932875000000002 +&
          q*(0.011359722222222243 + (0.088638125 -&
          0.030000000000000006*q)*q)) + 0.15915666666666684*Log(q)
    else if ((q < eps0).and.(q > -eps0)) then
      f = 1.
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
      df = 0.20915666666666688/q + q*(-0.9386575000000001 + q*(1.0340791666666667 + (-0.3954475 + 0.05*q)*q))
    else if (q >= 0.) then
      df = 0.15915666666666684/q + q*(-0.43865750000000003 + q*(0.03407916666666673 + (0.3545525 - 0.15000000000000002*q)*q))
    else if ((q < eps0).and.(q > -eps0)) then
      df = 0.
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
      ddf = -0.9386575000000001 - 0.20915666666666688/q**2 +&
            q*(1.0340791666666667 + (-0.3954475 + 0.05*q)*q) +&
            q*(1.0340791666666667 + (-0.3954475 + 0.05*q)*q + (-0.3954475 + 0.1*q)*q)
    else if (q >= 0.) then
      ddf = -0.43865750000000003 - 0.15915666666666684/q**2 +&
            q*(0.03407916666666673 + (0.3545525 - 0.15000000000000002*q)*q) +&
            q*(0.03407916666666673 + (0.3545525 - 0.30000000000000004*q)*q + (0.3545525 - 0.15000000000000002*q)*q)
    else if ((q < eps0).and.(q > -eps0)) then
      ddf = -3.
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine kddf
end module
