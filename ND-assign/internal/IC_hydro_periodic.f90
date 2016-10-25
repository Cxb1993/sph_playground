! subroutine periodic_ic(Xa, Xb, n, bn, rho, sk, pos, mas, vel, acc, den, slen, u)
!   integer, intent(in) :: n, bn
!   real, intent(in)    :: Xa, Xb, rho, sk
!   real, intent(out)   :: pos(n), vel(n), acc(n), mas(n), den(n), slen(n), u(n)
!   
!   integer :: i, nr
!   real :: step
!
!   nr = n - 2 * bn
!   step = (Xb - Xa) / nr
!
!   do i = 1, n
!     pos(i) = step * (i - bn)
!     mas(i) = step * rho
!     vel(i) = 0.0001 * sin(2.*pi*pos(i))
!     den(i) = rho
!     slen(i) = step * sk
!     acc(i) = 0.
!     u(i) = 1.
!   end do
!   ! 1   2   3   4   5   6   7   8   9   10  ......... 101 102 103 104 105 106 107 108 109 110
!   !                     106 107 108 109 110 ......... 1   2   3   4   5
!   do i = 1, bn
!     mas(i) = mas(nr + i)
!     vel(i) = vel(nr + i)
!     den(i) = den(nr + i)
!     slen(i) = slen(nr + i)
!
!     mas(nr + bn + i) = mas(bn + i)
!     vel(nr + bn + i) = vel(bn + i)
!     den(nr + bn + i) = den(bn + i)
!     slen(nr + bn + i) = slen(bn + i)
!   end do
! end subroutine periodic_ic
