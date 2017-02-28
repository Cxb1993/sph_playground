module err_calc
  use const
  use omp_lib
  use BC
  use kernel, only: get_dim,&
                    get_tasktype

  implicit none

  public :: err_init, err_T0sxsyet, err_infplate, err_sinxet,&
            err_diff_laplace, err_diff_graddiv, setStepsize

  private
  save
  real, allocatable :: tsin(:)
  real    :: period
  integer :: stepsize = 1

contains
  subroutine setStepsize(i)
    integer, intent(in) :: i
    stepsize = i
  end subroutine setStepsize

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

    print*, 'Not ready to NBS will divide to random number'

    !$OMP PARALLEL
    !$OMP DO PRIVATE(exact)
    do i=1,n
      exact = sin(pi*(pos(1,i)+1.)/2.) * sin(pi*(pos(2,i)+1.)/2.) * exp(-2.*(pi/2.)**2 * 0.1 * t)
      err(i) = abs(exact - num(i))
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine err_T0sxsyet

  subroutine err_sinxet(ptype, num, t, err, count)
    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(in)    :: num(:)
    real, allocatable, intent(inout) :: err(:)
    integer, intent(inout)           :: count
    real, intent(in)                 :: t

    integer             :: i, n
    real                :: exact

    print*, 'Not ready to NBS will divide to random number'

    n = size(ptype)
    count = 0
    err(1:n) = 0.
    !$omp parallel do default(none) &
    !$omp shared(n,tsin,num,err,t,ptype) &
    !$omp private(exact, i)&
    !$omp reduction(+:count)
    do i=1,n, stepsize
      if (ptype(i) /= 0) then
        exact = tsin(i) * exp(-pi**2 * t)
        err(i) = (exact - num(i))**2
        count = count + 1
      end if
    end do
    !$omp end parallel do
  end subroutine err_sinxet

  subroutine err_diff_laplace(ptype, x, num, err, count)
    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(in)    :: x(:,:), num(:,:)
    real, allocatable, intent(inout) :: err(:)
    integer, intent(inout)           :: count

    integer             :: i, n, dim
    real                :: exact(1:3)

    call get_dim(dim)
    n =size(ptype)
    count = 0
    err(:) = 0.
    exact(:) = 0.
    !$omp parallel do default(none) &
    !$omp shared(n,x,num,err,dim,period,ptype) &
    !$omp private(exact, i) &
    !$omp reduction(+:count)
    do i=1,n,stepsize
      exact(:) = 0.
      if (ptype(i) /= 0) then
        exact(1)  = -sin(period*x(1,i))
        if (dim > 1) then
          exact(2) = -sin(period*x(2,i))
        end if
        if (dim == 3) then
          exact(3) = -sin(period*x(3,i))
        end if
        err(i) = dot_product(exact(:) - num(:,i),exact(:) - num(:,i))
        count = count + 1
      end if
    end do
    !$omp end parallel do
  end subroutine err_diff_laplace

  subroutine err_diff_graddiv(ptype, x, num, err, count)
    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(in)    :: x(:,:), num(:,:)
    real, allocatable, intent(inout) :: err(:)
    integer, intent(out)             :: count

    integer             :: n, i, dim
    real                :: exact(1:3), xk(3)

    call get_dim(dim)
    n = size(ptype)
    count = 0
    err(:) = 0.
    !$omp parallel do default(none) &
    !$omp shared(n,ptype, x,num,err,dim,period) &
    !$omp private(exact, i,xk) &
    !$omp reduction(+:count)
    do i=1,n,stepsize
      ! print*,i, size(exact), size(x,dim=2),size(x,dim=1),size(num,dim=2), size(err), size(ptype)
      if (ptype(i) /= 0) then
        exact(:) = 0.
        if (dim == 1) then
          exact(1) = -(x(1,i))*Cos(x(1,i)) - 2*Sin(x(1,i))
        end if
        if (dim == 2) then
          exact(1) = -x(2,i)*Cos(x(1,i)) - Sin(x(2,i))
          exact(2) = -x(1,i)*Cos(x(2,i)) - Sin(x(1,i))
        end if
        if (dim == 3) then
          exact(1) = -(x(3,i)*Cos(x(1,i))) - Sin(x(2,i))
          exact(2) = -(x(1,i)*Cos(x(2,i))) - Sin(x(3,i))
          exact(3) = -(x(2,i)*Cos(x(3,i))) - Sin(x(1,i))
        end if
        err(i) = dot_product(exact(:)-num(:,i),exact(:)-num(:,i))
        count = count + 1
      end if
    end do
    !$omp end parallel do
  end subroutine err_diff_graddiv

  subroutine err_infplate(n, pos, num, t, err)
    integer, intent(in) :: n
    real, intent(in)    :: pos(3,n), num(n), t
    real, intent(inout) :: err

    integer :: i
    real :: tl, tr, tc, al, ar, ttmp, exact, xm, kl, kr, rhol, rhor, cvl, cvr

    print*, 'Not ready to NBS will divide to random number'

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
