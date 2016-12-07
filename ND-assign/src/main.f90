program main
  use BC
  use IC
  use kernel,   only:get_tasktype
  use iterator, only:iterate
  use printer
  use err_calc
  use args,     only:fillargs
  use bias,     only:get_chi
  use circuit1, only:c1_init

  implicit none

  real                :: sk = 1.2
  real                :: cv = 1.
  real                :: gamma = 1.4
  integer             :: n, dim
  character (len=40)  :: itype, errfname, ktype

  real                              :: dt, t, dtout, ltout, tfinish, error(11), npic, pspc1, pspc2
  integer                           :: iter, tt
  real, allocatable, dimension(:,:) :: p, v, a
  real, allocatable, dimension(:,:) :: pos, vel, acc, chi
  real, allocatable, dimension(:)   :: den, prs, mas, iu, du, o, c, h, dh, &
                                       cf, dcf, kcf, tdu, tdh, tcf, err

  print *, '#######################'
  print *, '#####'

  call fillargs(dim, pspc1, pspc2,&
                itype, ktype, errfname, dtout, npic, tfinish)

  sk = merge(1.2, -tfinish, tfinish > 0)
  tfinish = merge(tfinish, -1., tfinish > 0)

  call setup(itype, ktype, dim, n, sk, gamma, cv, &
                pspc1, pspc2, pos, vel, acc, &
                mas, den, h, prs, iu, du, cf, kcf, dcf)

  print *, '#####'
  print *, '#######################'

  call get_tasktype(tt)

  error(1) = pspc1
  error(2) = n

  t = 0.
  dt = 0.
  ltout = 0.
  p = pos
  v = vel
  a = acc

  allocate(err(1:n))
  allocate(chi(3,n))
  allocate(c(n))
  c(:) = cv
  allocate(tdu(n))
  allocate(dh(n))
  dh(:) = 0
  allocate(tdh(n))
  allocate(tcf(n))

  read *

  call err_init(n, pos)
  call c1_init(n)

  do while (t <= tfinish)
    print *, t
    select case(tt)
    case(1)
      ! 'hydroshock'
      dt = .3 * minval(h) / maxval(c)
    case(2)
      ! 'infslb'
      dt = .144 * minval(den) * minval(c) * minval(h) ** 2 / maxval(kcf)
    case(3)
      ! 'hc-sinx'
      dt = .144 * minval(den) * minval(c) * minval(h) ** 2 / maxval(kcf)
      call err_sinxet(n, cf, t, err)
    end select

    if (t >= ltout) then
      print *, iter, t, sqrt(sum(err)/n)
      call print_output(n, t, pos, vel, acc, mas, den, h, prs, iu, cf, sqrt(err))
      ltout = ltout + dtout
    end if

    p(:,:) = pos(:,:)
    v(:,:) = vel(:,:)
    a(:,:) = acc(:,:)
    tdu(:) = du(:)
    tdh(:) = dh(:)
    tcf(:) = dcf(:)

    pos(:,:) = p(:,:) + dt * v(:,:) + 0.5 * dt * dt * a(:,:)
    vel(:,:) = v(:,:) + dt * a(:,:)
    iu(:)   = iu(:) + dt * du(:)
    h(:)     = h(:)   + dt *  dh(:)
    ! cf(:)    = cf(:)  + dt * dcf(:)

    call iterate(n, sk, gamma, &
                pos, vel, acc, &
                mas, den, h, dh, o, prs, c, iu, du, &
                cf, dcf, kcf)

    vel(:,:) = vel(:,:) + 0.5 * dt * (acc(:,:) - a(:,:))
    iu(:)   = iu(:)   + 0.5 * dt * (du(:) - tdu(:))
    h(:)     = h(:)     + 0.5 * dt * (dh(:) - tdh(:))

    select case(tt)
    case(1)
      ! 'hydroshock'
      cf(:) = cf(:)     + 0.5 * dt * (dcf(:) - tcf(:))
    case(2, 3)
      ! 'infslb', 'hc-sinx'
      cf(:) = iu(:) / cv
    end select

    t = t + dt
    iter = iter + 1

    if ((t-dt/2<tfinish*1/3).and.(tfinish*1/3<t+dt/2)) then
      error(3) = sqrt(sum(err)/n)
      call get_chi(n, mas, den, pos, h, chi)
      call periodic1(chi, 1)
      error(4) = sum(chi)/n/dim
      error(5) = t
    else if ((t-dt/2<tfinish*2/3).and.(tfinish*2/3<t+dt/2)) then
      error(6) = sqrt(sum(err)/n)
      call get_chi(n, mas, den, pos, h, chi)
      error(7) = sum(chi)/n/dim
      error(8) = t
    end if
  end do

  if (t == .0) then
    select case(tt)
    case(1)
      ! 'hydroshock'
    case(2)
      ! 'infslb'
    case(3)
      ! 'hc-sinx'
      call err_sinxet(n, cf, t, err)
    end select
    call get_chi(n, mas, den, pos, h, chi)
    call periodic3(chi, 00, 0)
    call print_output(n, t, pos, chi, acc, mas, den, h, prs, iu, cf, sqrt(err))
    error(3) = sqrt(sum(err)/n)
    error(4) = sum(chi)/n/dim
    error(5) = sk

    call print_appendline(5, error, errfname)
  else
    error(9) = sqrt(sum(err)/n)
    call get_chi(n, mas, den, pos, h, chi)
    error(10) = sum(chi)/n/dim
    error(11) = t
    call print_appendline(11, error, errfname)
  end if
end program main
