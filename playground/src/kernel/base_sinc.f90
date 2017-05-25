module sinc
  use const

  implicit none

  public :: kf, kdf, kddf, wCv, krad, kernelname, setdimbase

  private

    real :: knorm(3) = (/ 9.45261263832083, 7.29241786955192, 5.76692759942686 /)
    real :: krad = 2., wCv
    integer :: dim
    character (len=10) :: kernelname='sinc'

 contains
   subroutine setdimbase(d)
     integer, intent(in) :: d
     dim = d
     wCv = knorm(dim)
   end subroutine

  pure subroutine kf(r, h, f)
    real, intent(in)  :: r, h
    real, intent(out) :: f
    real              :: q

    q = r / h
    if (q >= 2.) then
      f  = 0.
    else if (q > 0.) then
      f = 4*sin(pi*q/2)**4/(pi**5*q**4)
    else if ((q < epsilon(1)).and.(q > -epsilon(1))) then
      f = 0.0795775
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
    if (q >= 2.) then
      df = 0.
    else if (q > 0.) then
      df = 8*(pi*q*cos(pi*q/2) - 2*sin(pi*q/2))*sin(pi*q/2)**3/(pi**5*q**5)
    else if ((q < epsilon(1)).and.(q > -epsilon(1))) then
      f = 0.
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
    if (q >= 2.) then
      ddf = 0.
    else if (q >= 0.) then
      ddf = 4*(2*pi**2*q**2*cos(pi*q) + pi**2*q**2 - 8*pi*q*sin(pi*q) - 10*cos(pi*q) + 10)*sin(pi*q/2)**2/(pi**5*q**6)
    else
      if (isnan(q)) then
        error stop 'q is nan'
      else
        error stop 'q is negative'
      end if
    end if
  end subroutine kddf
end module sinc
