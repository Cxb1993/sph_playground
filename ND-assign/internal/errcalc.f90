module err_calc
  use omp_lib
  use BC

  implicit none

  public :: err_T0sxsyet, err_infplate, err_nonhplate

  private
  real, parameter     :: pi = 4.*atan(1.)

contains
  subroutine err_T0sxsyet(n, pos, num, t, err)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n,3), num(n), t
    real, intent(out)   :: err(n)

    integer             :: i
    real                :: exact

    !$OMP PARALLEL
    !$OMP DO PRIVATE(exact)
    do i=1,n
      exact = sin(pi*(pos(i,1)+1.)/2.) * sin(pi*(pos(i,2)+1.)/2.) * exp(-2.*(pi/2.)**2 * 0.1 * t)
      if (exact < 1E-10) then
        err(i) = 0.
      else
        err(i) = abs((exact - num(i)) / exact)
      end if
      ! print *, exact, num(i), err(i), err(i) * 100
      ! read *
      ! err(i) = abs(sin(pi*(pos(i,1)+1.)/2.) * exp(-(pi/2.)**2 * 0.1 * t) - num(i))
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    ! print *, ''
    ! print *, ''
    ! print *, maxval(err), sum(err)/size(err), sum(err)/size(err) * 100
    ! print *, ''
    ! read *

  end subroutine err_T0sxsyet

  subroutine err_infplate(n, pos, num, t, err)
    integer, intent(in) :: n
    real, intent(in) :: pos(n,3), num(n), t
    real, intent(inout) :: err(n)

    integer :: i
    real :: tl, tr, tc, al, ar, ttmp, exact, xm, kl, kr, rhol, rhor, cvl, cvr

    kl = 1000.
    kr = 1.
    rhol = 1.
    rhor = 1.
    cvl = 1.
    cvr = 1.
    tl = 1.
    tr = 2.
    xm = 0.

    al = kl / rhol / cvl
    ar = kr / rhor / cvr

    tc = (tr - tl) * (kr / sqrt(ar)) / (kr / sqrt(ar) + kl / sqrt(al))
    !$OMP PARALLEL
    !$OMP DO PRIVATE(exact, ttmp)
    do i=1,n
      if (pos(i,1) < xm) then
        ttmp = erfc((xm-pos(i,1))/(2 * sqrt(al*t)))
      else
        ttmp = 1 + (kl/kr)*sqrt(ar/al)*erf((pos(i,1)-xm)/(2 * sqrt(ar*t)))
      end if
      exact = ttmp * tc + tl
      ! err(i) = ttmp * tc + tl
      if (abs(exact - num(i)) < 1E-10) then
        err(i) = 0.
      else
        err(i) = abs((exact - num(i)) / exact)
      end if
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    ! trn = (1 + (kl/kr)*sqrt(ar/al)*erf((xr-xm)/sqrt(ar*time))) * tc + tln
    ! tln = erfc((xm-xl)/(al*time)) * tc + tln
    return
  end subroutine err_infplate

  subroutine err_nonhplate(n, pos, num, t, err)
    integer, intent(in) :: n
    real, intent(in) :: pos(n,3), num(n), t
    real, intent(inout) :: err(n)

    integer :: i, j
    real :: tl, tr, tc, al, ar, ttmp, exact, xm, kl, kr, rhol, rhor, cvl, cvr, f, x, li

    kl = 10.
    kr = 1.
    rhol = 1.
    rhor = 1.
    cvl = 1.
    cvr = 1.
    tl = 4.
    tr = 1.
    xm = 0.

    al = kl / rhol / cvl
    ar = kr / rhor / cvr
    tc = (tl - tr) * (kr / sqrt(ar)) / (kr / sqrt(ar) + kl / sqrt(al))

    do i=1,n
      ttmp = 0.
      if (pos(i,1) < xm) then
        x = pos(i,1)/(xm-pos(1,1))
        f = al*t/(xm-pos(1,1))**2
        do j=1,100
          li = j*pi
          ttmp = ttmp - 2./pi * sin(li * x) * exp(-li**2*f) / j
        end do
        ttmp = 1 + (kl/kr)*sqrt(ar/al)*(ttmp - x)
        exact = ttmp*tc + tr
      else
        x = pos(i,1)/(pos(n,1)-xm)
        f = ar*t/(pos(n,1)-xm)**2
        do j=1,100
          li = j*pi
          ttmp = ttmp - 2./pi * sin(li * x) * exp(-li**2*f) / j
        end do
        ttmp = 1 - x + ttmp
        exact = ttmp * tc + tr
      end if
      ! print *, exact
      err(i) = exact
    end do
    ! read *
  end subroutine err_nonhplate
end module err_calc
