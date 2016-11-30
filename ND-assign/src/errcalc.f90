module err_calc
  use omp_lib
  use BC
  use kernel

  implicit none

  public :: err_init, err_T0sxsyet, err_infplate, err_sinxet

  private
  save
  real, parameter   :: pi = 4.*atan(1.)
  real, allocatable :: tsin(:)

contains
  subroutine err_init(n, pos)
    integer, intent(in) :: n
    real, intent(in)    :: pos(3,n)
    integer             :: tt, i

    allocate(tsin(n))

    call get_tasktype(tt)

    select case(tt)
    case(1)
      ! 'hydroshock'
    case(2)
      ! 'infslb'
    case(3)
      ! 'hc-sinx'
      do i=1,n
        tsin(i) = sin(pi * (pos(1,i) + 1.))
      end do
    end select


  end subroutine err_init

  subroutine err_T0sxsyet(n, pos, num, t, err)
    integer, intent(in) :: n
    real, intent(in)    :: pos(3,n), num(n), t
    real, intent(out)   :: err(n)

    integer             :: i
    real                :: exact

    !$OMP PARALLEL
    !$OMP DO PRIVATE(exact)
    do i=1,n
      exact = sin(pi*(pos(1,i)+1.)/2.) * sin(pi*(pos(2,i)+1.)/2.) * exp(-2.*(pi/2.)**2 * 0.1 * t)
      err(i) = abs(exact - num(i))
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine err_T0sxsyet

  subroutine err_sinxet(n, num, t, err)
    integer, intent(in) :: n
    real, intent(in)    :: num(n), t
    real, intent(out)   :: err(n)

    integer             :: i
    real                :: exact

    !$omp parallel do default(none) &
    !$omp shared(n,pos,num,err,t) &
    !$omp private(exact, i)
    do i=1,n
      exact = tsin(i) * exp(-pi**2 * t)
      ! exact = sin(pi * (pos(1,i) + 1.)) * exp(-pi**2 * t)
      err(i) = (exact - num(i))**2
    end do
    !$omp end parallel do
  end subroutine err_sinxet

  subroutine err_infplate(n, pos, num, t, err)
    integer, intent(in) :: n
    real, intent(in)    :: pos(3,n), num(n), t
    real, intent(inout) :: err

    integer :: i
    real :: tl, tr, tc, al, ar, ttmp, exact, xm, kl, kr, rhol, rhor, cvl, cvr

    kl = 1.
    kr = 1.
    rhol = 1.
    rhor = 1.
    cvl = 1.
    cvr = 1.
    tl = 0.
    tr = 1.
    xm = 0.

    al = kl / rhol / cvl
    ar = kr / rhor / cvr

    tc = (tr - tl) * (kr / sqrt(ar)) / (kr / sqrt(ar) + kl / sqrt(al))

    err = 0
    !$omp parallel do default(none) &
    !$omp shared(pos,n,xm,al,ar,kl,kr,tc,tl,num,t) &
    !$omp private(exact, ttmp, i) &
    !$omp reduction(+:err)
    do i=1,n
      if (pos(1,i) < xm) then
        ttmp = erfc((xm-pos(1,i))/(2 * sqrt(al*t)))
      else
        ttmp = 1 + (kl/kr)*sqrt(ar/al)*erf((pos(1,i)-xm)/(2 * sqrt(ar*t)))
      end if
      exact = ttmp * tc + tl
      err = err + (num(i) - exact)**2
    end do
    !$omp end parallel do
    err = sqrt(err/n)
    return
  end subroutine err_infplate
end module err_calc
