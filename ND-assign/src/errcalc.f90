module err_calc
  use const
  use omp_lib
  use BC
  use kernel

  implicit none

  public :: err_init, err_T0sxsyet, err_infplate, err_sinxet, err_diff_laplace, err_diff_graddiv

  private
  save
  real, allocatable :: tsin(:)
  real :: period

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
        ! print *, pos(:,i)
        tsin(i) = sin(pi * (pos(1,i) + 1.))
      end do
    case(5)
      ! 'diff-laplace'
      period = 1.
    case(6)
      ! 'diff-graddiv
      period = 1.
    case default
      print *, 'Task type was not sen in error init errcalc.f90'
      stop
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

    ! print *, maxval(num)
    !$omp parallel do default(none) &
    !$omp shared(n,tsin,num,err,t) &
    !$omp private(exact, i)
    do i=1,n
      exact = tsin(i) * exp(-pi**2 * t)
      ! exact = sin(pi * (pos(1,i) + 1.)) * exp(-pi**2 * t)
      ! print *, maxval(num), minval(num), num(i), n
      err(i) = (exact - num(i))**2
    end do
    !$omp end parallel do
  end subroutine err_sinxet

  subroutine err_diff_laplace(n, x, num, dim, err)
    integer, intent(in) :: n, dim
    real, intent(in)    :: x(3,n), num(3,n)
    real, intent(out)   :: err(n)

    integer             :: i
    real                :: exact(3)

    !$omp parallel do default(none) &
    !$omp shared(n,x,num,err,dim,period) &
    !$omp private(exact, i)
    do i=1,n
      exact(:) = sin(period*x(:,i))
      ! if (dim > 1) then
      !   exact = exact * sin(period*x(2,i))
      !   if (dim == 3) then
      !     exact = exact * sin(period*x(3,i))
      !   end if
      ! end if
      ! exact(:) = -dim*period**2*exact(:)
      exact(:) = -period**2*exact(:)
      err(i) =dot_product(exact - num(:,i),exact - num(:,i))
      ! print *, i, exact, num(1,i), err(i)
    end do
    !$omp end parallel do
  end subroutine err_diff_laplace

  subroutine err_diff_graddiv(n, x, num, dim, err)
    integer, intent(in) :: n, dim
    real, intent(in)    :: x(3,n), num(3,n)
    real, intent(out)   :: err(n)

    integer             :: i
    real                :: exact(3)

    !$omp parallel do default(none) &
    !$omp shared(n,x,num,err,dim,period) &
    !$omp private(exact, i)
    do i=1,n
      exact(:) = sin(period*x(:,i))
      ! if (dim > 1) then
      !   exact = exact * sin(period*x(2,i))
      !   if (dim == 3) then
      !     exact = exact * sin(period*x(3,i))
      !   end if
      ! end if
      exact(:) = -exact(:)
      err(i) = dot_product(exact(:)-num(:,i),exact(:)-num(:,i))
      ! print *, i
      ! print *, exact(:)
      ! print *, num(:,i)
      ! print *, err(i)
      ! read *
    end do
    !$omp end parallel do
  end subroutine err_diff_graddiv

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
