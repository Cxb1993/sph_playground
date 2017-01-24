program main
  use BC
  use IC,       only:setupIC
  use kernel,   only:get_tasktype
  use iterator, only:iterate
  use printer
  use err_calc, only:err_init,&
                     err_diff_laplace,&
                     err_diff_graddiv,&
                     err_sinxet
  use args,     only:fillargs
  use bias,     only:get_chi
  use circuit1, only:c1_init

  implicit none

  real                :: sk = 1.2
  real                :: cv = 1.
  integer             :: n, dim, i
  character (len=40)  :: itype, errfname, ktype

  real                              :: dt, t, dtout, ltout, tfinish, error(11), npic,&
                                       pspc1, pspc2, gamma
  integer                           :: iter, tt
  real, allocatable, dimension(:,:) :: p, v, a
  real, allocatable, dimension(:,:) :: pos, vel, acc, chi
  real, allocatable, dimension(:)   :: den, prs, mas, iu, du, om, c, h, dh, &
                                       cf, dcf, kcf, tdu, tdh, tcf, err

  print *, '#######################'
  print *, '#####'

  call fillargs(dim, pspc1, pspc2,&
                itype, ktype, errfname, dtout, npic, tfinish)

  sk = merge(sk, -tfinish, tfinish > 0)
  tfinish = merge(tfinish, -1., tfinish > 0)

  call setupIC(n, sk, gamma, cv, pspc1, pspc2, pos, vel, acc, &
                mas, den, h, prs, iu, du, cf, kcf, dcf)
  ! print *, 0, -2
  ! print *, maxval(abs(du))

  print *, '#####'
  print *, '#######################'

  call get_tasktype(tt)

  error(1) = pspc1
  error(2) = n

  t = 0.
  dt = 0.
  ltout = 0.
  iter = 0.

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
  allocate(om(n))

  read *

  call err_init(n, pos)
  call c1_init(n)

 print *, tfinish
  do while (t <= tfinish)
    ! print *, 0, -1
    ! print *, '--0'
    ! print *, t
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
    case(4)
      ! 'photoevaporation' 'pheva'
      dt = .3e-3 * minval(h)**2 / maxval(c)**2 / maxval(kcf) / maxval(cf)
    case (5,6)
      ! 'diff-laplass'      ! 'diff-graddiv'
      dt = 0.
      tfinish = -1.
    case default
      print *, 'Task type time increment was not defined'
      stop
    end select
    ! print *, 0, 0
    if (t >= ltout) then
      print *, iter, t, dt
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
    iu(:)    = iu(:) + dt * du(:)
    h(:)     = h(:)   + dt *  dh(:)
    if (tt == 4) then
      cf(:)    = cf(:)  + dt * dcf(:)
    end if
    ! print *, 0, 1
    ! print *, maxval(abs(du))
    call iterate(n, sk, gamma, &
                pos, vel, acc, &
                mas, den, h, dh, om, prs, c, iu, du, &
                cf, dcf, kcf)
    ! print *, maxval(abs(du))
    ! print *, 0, 2
    ! print *, maxval(cf)

    vel(:,:) = vel(:,:) + 0.5 * dt * (acc(:,:) - a(:,:))
    iu(:)   = iu(:)   + 0.5 * dt * (du(:) - tdu(:))
    h(:)     = h(:)     + 0.5 * dt * (dh(:) - tdh(:))

    select case(tt)
    case(1, 4, 5, 6)
      ! 'hydroshock' ! 'diff-graddiv'
      cf(:) = cf(:)     + 0.5 * dt * (dcf(:) - tcf(:))
    case(2, 3)
      ! 'infslb', 'hc-sinx'
      cf(:) = iu(:) / cv
    case default
      print *, 'Task type was not sen in couples field integration'
      stop
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
    ! print *, '--1'
    ! print *, 0, 4
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
      call get_chi(n, mas, den, pos, h, chi)
      call periodic3(chi, 00, dim)
      call print_output(n, t, pos, chi, acc, mas, den, h, prs, iu, cf, sqrt(err))
      error(4) = sum(chi)/n/dim
    case(5)
      ! 'diff-laplace'
      call err_diff_laplace(n, pos, acc, dim, err)
      call get_chi(n, mas, den, pos, h, chi)
      call periodic3(chi, 00, dim)
      call print_output(n, t, pos, chi, acc, mas, den, h, prs, iu, cf, sqrt(err))
      error(4) = sum(chi)/n/dim
    case(6)
      ! 'diff-graddiv'
      call err_diff_graddiv(n, pos, acc, dim, err)
      call get_chi(n, mas, den, pos, h, chi)
      call periodic3(chi, 00, dim)
      call print_output(n, t, pos, chi, acc, mas, den, h, prs, iu, cf, sqrt(err))
      error(4) = sum(chi)/n/dim
    case default
      print *, 'Task type was not sen in error evaluation main.f90'
      stop
    end select
    error(3) = sqrt(sum(err)/n)
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
