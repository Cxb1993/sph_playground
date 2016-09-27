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

  subroutine shock_ic(dim, nx, nb, sk, g, pos, vel, acc, mas, den, sln, p, u, n)
    integer, intent(in)  :: dim, nx, nb
    real, intent(in)     :: sk, g
    real, intent(out)    :: pos(dim,nx), mas(dim,nx), vel(dim,nx), den(dim,nx), sln(dim,nx), acc(dim,nx)
    real, intent(out)    :: p(dim,nx), u(dim,nx)
    integer, intent(out) :: n(dim)
    integer              :: d, i
    real                 :: spatVarBrdrs(dim,2), parSpacing(2), shockPressure(2), shockDensity(2)
    real                 :: x, svb1, svb2, p1, p2, r1, r2, l1, l2

    spatVarBrdrs(1,1) = -0.5
    spatVarBrdrs(1,2) = 0.5
    spatVarBrdrs(2,1) = 0.
    spatVarBrdrs(2,2) = 0.1

    parSpacing(1) = 0.001
    parSpacing(2) = 0.008

    shockPressure(1) = 1.
    shockPressure(2) = 0.1

    shockDensity(1) = 1.
    shockDensity(2) = 0.125

    do d = 1,dim
      svb1 = spatVarBrdrs(d,1)
      svb2 = spatVarBrdrs(d,2)
      p1 = shockPressure(1)
      p2 = shockPressure(2)
      r1 = shockDensity(1)
      r2 = shockDensity(2)
      l1 = parSpacing(1)
      l2 = parSpacing(2)

      x = svb1
      i = 1
      do while ((x >= svb1).and.(x < svb2))
        print  *, x
        if (x<0) then
          pos(d,i) = x
          mas(d,i) = l1 * r1
          vel(d,i) = 0.
          den(d,i) = r1
          sln(d,i) = l1 * sk
          acc(d,i) = 0.
          p(d,i) = p1
          u(d,i) = p1 / (g - 1) / r1
          x = x + l1
        else
          pos(d,i) = x
          mas(d,i) = l2 * r2
          vel(d,i) = 0.
          den(d,i) = r2
          sln(d,i) = l2 * sk
          acc(d,i) = 0.
          p(d,i) = p2
          u(d,i) = p2 / (g - 1) / r2
          x = x + l2
        end if
      end do
    end do    !
    ! do i = 1, nb+1
    !   pos(na+i) = stepb * (i-1)
    !   mas(na+i) = stepb * rhob
    !   vel(na+i) = 0.
    !   den(na+i) = rhob
    !   slen(na+i) = stepb * sk
    !   acc(na+i) = 0.
    !   p(na+i) = pb
    !   u(na+i) = pb / (gamma - 1) / den(na+i)
    ! end do
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
