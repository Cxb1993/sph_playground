module setup
  implicit none

contains

  subroutine periodic(fromX, toX, n, bn, rho, sk, pos, mas, vel, den, slen)
    integer :: n, bn
    real, intent(in)    :: fromX, toX, rho, sk
    real, intent(out)   :: pos(n), mas(n), vel(n), den(n), slen(n)
    real, parameter :: pi = 4.*atan(1.)
    integer :: i
    real :: step

    step = (toX - fromX) / (n - 2 * bn - 1)
    do i = 1, n
      pos(i) = step * (i - bn - 1)
      mas(i) = (1.0 * rho)/(n - 2 * bn)
      vel(i) = sin(2.*pi*pos(i)) * (0.0001)
      den(i) = rho
      slen(i) = step * sk
    end do
    do i = 1, bn
      mas(bn - i + 1) = mas(n - bn - i)
      vel(bn - i + 1) = vel(n - bn - i)
      den(bn - i + 1) = den(n - bn - i)
      slen(bn - i + 1) = slen(n - bn - i)

      mas(n - bn + i) = mas(bn + i + 1)
      vel(n - bn + i) = vel(bn + i + 1)
      den(n - bn + i) = den(bn + i + 1)
      slen(n - bn + i) = slen(bn + i + 1)
    end do
  end subroutine periodic
end module setup
