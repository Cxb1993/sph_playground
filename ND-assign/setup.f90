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

  subroutine shock_ic(dim, nx, n, sk, g, pos, vel, acc, mas, den, sln, prs, uie)
    integer, intent(in)  :: dim, nx
    real, intent(in)     :: sk, g
    real, intent(out)    :: pos(3,nx), mas(nx), vel(nx), den(nx), sln(nx), acc(nx), prs(nx), uie(nx)
    integer, intent(out) :: n
    real                 :: spatVarBrdrs11, spatVarBrdrs12, spatVarBrdrs21, spatVarBrdrs22, spatVarBrdrs31, spatVarBrdrs32
    real                 :: parSpacing1, parSpacing2, shockPressure1, shockPressure2, shockDensity1, shockDensity2
    real                 :: x, y, z

    spatVarBrdrs11 = -0.5
    spatVarBrdrs12 = 0.5
    spatVarBrdrs21 = 0.
    spatVarBrdrs22 = 0.
    spatVarBrdrs31 = 0.
    spatVarBrdrs32 = 0.
    if (dim.gt.1) then
      spatVarBrdrs21 = 0.
      spatVarBrdrs22 = 0.1
    end if
    if (dim.eq.3) then
      spatVarBrdrs31 = 0.
      spatVarBrdrs32 = 0.1
    end if

    parSpacing1 = 0.001
    parSpacing2 = 0.008

    shockPressure1 = 1.
    shockPressure2 = 0.1

    shockDensity1 = 1.
    shockDensity2 = 0.125

    x = spatVarBrdrs11
    n = 1
    do while ((x >= spatVarBrdrs11).and.(x <= spatVarBrdrs12))
      y = spatVarBrdrs21
      do while ((y >= spatVarBrdrs21).and.(y <= spatVarBrdrs22))
        z = spatVarBrdrs31
        do while ((z >= spatVarBrdrs31).and.(z <= spatVarBrdrs32))
          pos(1,n) = x
          pos(2,n) = y
          pos(3,n) = z
          if (x<0) then
            mas(n) = (parSpacing1**3) * shockDensity1
            vel(n) = 0.
            den(n) = shockDensity1
            sln(n) = (parSpacing1**3) * sk
            acc(n) = 0.
            prs(n) = shockPressure1
            uie(n) = shockPressure1 / (g - 1) / shockDensity1
          else
            mas(n) = (parSpacing2**3) * shockDensity2
            vel(n) = 0.
            den(n) = shockDensity2
            sln(n) = (parSpacing2**3) * sk
            acc(n) = 0.
            prs(n) = shockPressure2
            uie(n) = shockPressure2 / (g - 1) / shockDensity2
          end if
          if (x<0) then
            z = z + parSpacing1
          else
            z = z + parSpacing2
          end if
          n = n + 1
        end do
        if (x<0) then
          y = y + parSpacing1
        else
          y = y + parSpacing2
        end if
      end do
      if (x<0) then
        x = x + parSpacing1
      else
        x = x + parSpacing2
      end if
    end do
    print *, 'Particles placed : ', n
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
