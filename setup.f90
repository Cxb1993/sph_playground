module setup
  implicit none

  public :: periodic_ic, set_periodic, shock_ic, set_fixed

  private

contains

  subroutine periodic_ic(Xa, Xb, n, bn, rho, sk, pos, mas, vel, acc, den, slen, u)
    integer, intent(in) :: n, bn
    real, intent(in)    :: Xa, Xb, rho, sk
    real, intent(out)   :: pos(n), mas(n), vel(n), den(n), slen(n), acc(n), u(n)
    real, parameter     :: pi = 4.*atan(1.)
    integer :: i, nr
    real :: step

    nr = n - 2 * bn
    step = (Xb - Xa) / nr

    do i = 1, n
      pos(i) = step * (i - bn)
      mas(i) = step * rho
      vel(i) = 0.0001 * sin(2.*pi*pos(i))
      den(i) = rho
      slen(i) = step * sk
      acc(i) = 0.
      u(i) = 1.
    end do
! 1   2   3   4   5   6   7   8   9   10  ......... 101 102 103 104 105 106 107 108 109 110
!                     106 107 108 109 110 ......... 1   2   3   4   5
    do i = 1, bn
      mas(i) = mas(nr + i)
      vel(i) = vel(nr + i)
      den(i) = den(nr + i)
      slen(i) = slen(nr + i)

      mas(nr + bn + i) = mas(bn + i)
      vel(nr + bn + i) = vel(bn + i)
      den(nr + bn + i) = den(bn + i)
      slen(nr + bn + i) = slen(bn + i)
    end do
  end subroutine periodic_ic

  subroutine set_periodic(n, bn, A)
    integer, intent(in) :: n, bn
    real, intent(out)   :: A(n)
    integer             :: i, nr

    nr = n - 2 * bn
    do i = 1, bn
      A(i) = A(nr + i)
      A(nr + bn + i) = A(bn + i)
    end do
  end subroutine set_periodic

  subroutine shock_ic(Xa, Xb, na, nb, rhoa, rhob, pa, pb, sk, pos, mas, vel, acc, den, slen, p, u, gamma)
    integer, intent(in) :: na, nb
    real, intent(in)    :: Xa, Xb, rhoa, rhob, sk, pa, pb, gamma
    real, intent(out)   :: pos(na+nb+1), mas(na+nb+1), vel(na+nb+1), den(na+nb+1), slen(na+nb+1), acc(na+nb+1)
    real, intent(out)   :: p(na+nb+1), u(na+nb+1)
    integer             :: i
    real                :: stepa, stepb

    stepa = Xa / na
    stepb = Xb / nb

    do i = 1, na
      pos(i) = -(Xa - stepa * (i-1))
      mas(i) = stepa * rhoa
      vel(i) = 0.
      den(i) = rhoa
      slen(i) = stepa * sk
      acc(i) = 0.
      p(i) = pa
      u(i) = pa / (gamma - 1) / den(i)
    end do

    do i = 1, nb+1
      pos(na+i) = stepb * (i-1)
      mas(na+i) = stepb * rhob
      vel(na+i) = 0.
      den(na+i) = rhob
      slen(na+i) = stepb * sk
      acc(na+i) = 0.
      p(na+i) = pb
      u(na+i) = pb / (gamma - 1) / den(na+i)
    end do
  end subroutine shock_ic

  subroutine set_fixed(n, bn, A)
    integer, intent(in) :: n, bn
    real, intent(out)   :: A(n)
    integer             :: i

    do i = 1, bn
      A(i) = 0.
      A(n-i+1) = 0.
    end do
  end subroutine set_fixed

end module setup
