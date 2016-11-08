program main
  use BC
  use IC
  use internal
  use printer
  use kernel
  use err_calc

  implicit none

  character (len=40) :: itype, snpfname, errfname
  character (len=1)  :: arg

  real                :: sk = 1.2
  integer, parameter  :: nmax = 400000
  real                :: cv = 1.
  real                :: gamma = 1.4
  integer             :: n, dim

  real, allocatable, dimension(:,:) :: position, velocity, acceleration
  real :: density(nmax), slength(nmax), pressure(nmax), mass(nmax), ienergy(nmax), dienergy(nmax), omega(nmax)
  real :: coupledfield(nmax), kcoupledfield(nmax), dcoupledfield(nmax)

  real                              :: dt, t, dtout, ltout, tfinish, error(5),&
                                       pspc1, pspc2, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2
  integer                           :: finish, iter, npic
  real, allocatable, dimension(:,:) :: p, v, a
  real, allocatable, dimension(:,:) :: pos, vel, acc
  real, allocatable, dimension(:)   :: den, prs, mas, ieu, diu, o, du, c, h, dh, tdh
  real, allocatable, dimension(:)   :: cf, tcf, dcf, kcf

  allocate(position(nmax,3))
  allocate(velocity(nmax,3))
  allocate(acceleration(nmax,3))

  print *, '####################'
  print *, '#'
  call get_command_argument(1, arg)
  read(arg(:), fmt="(i5)") dim
  print *, "#       dim:", dim

  call get_command_argument(2, itype)
  print *, "# task type:   ", itype

  pspc1 = 1./5
  pspc2 = pspc1
  print *, "# p.spacing:   ", pspc1, pspc2
  brdx1 = -10.
  brdx2 = 10.
  if (dim.gt.1) then
    brdy1 = - pspc1 * 4.5
    brdy2 = pspc2 * 4.5
  else
    brdy1 = 0.
    brdy2 = 0.
  end if
  if (dim.eq.3) then
    brdz1 = - pspc1 * 4.5
    brdz2 = pspc2 * 4.5
  else
    brdz1 = 0.
    brdz2 = 0.
  end if
  call setup(itype, dim, nmax, n, sk, gamma, cv, &
                pspc1, pspc2, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, &
                position, velocity, acceleration, &
                mass, density, slength, pressure, ienergy, coupledfield, kcoupledfield, dcoupledfield)

  error(1) = pspc1
  error(2) = n
  snpfname = 'n2wc-infslb-1-2D-snp'
  errfname = 'n2wc-infslb-1-2D-err'

  t = 0.
  dt = 0.
  ltout = 0.
  finish = 3000.
  tfinish = 5.
  npic = 200.
  dtout = tfinish / npic
  print *, '#  print dt:', dtout
  print *, '#'
  print *, '####################'

  read *

  pos = position(1:n,:)
  p   = pos(:,:)
  vel = velocity(1:n,:)
  v   = vel(:,:)
  acc = acceleration(1:n,:)
  a   = acc(:,:)
  deallocate(position)
  deallocate(velocity)
  deallocate(acceleration)
  den = density(1:n)
  h   = slength(1:n)
  prs = pressure(1:n)
  mas = mass(1:n)
  ieu = ienergy(1:n)
  diu = dienergy(1:n)
  du  = diu(:)
  o   = omega(1:n)
  allocate(c(1:n))
  c(:)= cv
  allocate(dh(1:n))
  dh  = 0
  tdh = dh
  cf  = coupledfield(1:n)
  tcf = coupledfield(1:n)
  allocate(dcf(1:n))
  dcf = dcoupledfield(1:n)
  kcf = kcoupledfield(1:n)

  call derivs(n, sk, gamma, &
              pos, vel, acc, &
              mas, den, h, dh, o, prs, c, ieu, diu, &
              cf, dcf, kcf)
  call print_output(n, 0., pos, vel, acc, mas, den, h, prs, ieu, cf)

  call plot_simple(2, error(1:2), snpfname)
  do while (t <= tfinish)
    select case(itype)
    case('hydroshock')
      dt = .3 * minval(h) / maxval(c)
    case('heatslab')
      dt = .144 * minval(den) * minval(c) * minval(h) ** 2 / maxval(kcf)
    end select
    ! print *, "dt: ", dt, minval(den), minval(c), minval(h), maxval(kcf)
    ! read *
    if (t >= ltout) then
      print *, iter, t
      call print_output(n, t, pos, vel, acc, mas, den, h, prs, ieu, cf)
      ltout = ltout + dtout
    end if
    p(:,:) = pos(:,:)
    v(:,:) = vel(:,:)
    a(:,:) = acc(:,:)
    du(:)  = diu(:)
    tdh(:) = dh(:)
    tcf(:) = dcf(:)

    pos(:,:) = p(:,:) + dt * v(:,:) + 0.5 * dt * dt * a(:,:)
    vel(:,:) = v(:,:) + dt * a(:,:)
    ieu(:)   = ieu(:) + dt * diu(:)
    h(:)     = h(:)   + dt *  dh(:)
    ! cf(:)    = cf(:)  + dt * dcf(:)

    call derivs(n, sk, gamma, &
                pos, vel, acc, &
                mas, den, h, dh, o, prs, c, ieu, diu, &
                cf, dcf, kcf)

    vel(:,:) = vel(:,:) + 0.5 * dt * (acc(:,:) - a(:,:))
    ieu(:)   = ieu(:)   + 0.5 * dt * (diu(:) - du(:))
    h(:)     = h(:)     + 0.5 * dt * (dh(:) - tdh(:))
    if (itype.eq.'heatslab') then
      cf(:) = ieu(:) / cv
    else
      cf(:) = cf(:)     + 0.5 * dt * (dcf(:) - tcf(:))
    end if

    if ((t-dt/2<tfinish*1/3).and.(tfinish*1/3<t+dt/2)) then
      call plot_simple(n, pos, snpfname)
      call plot_simple(n, cf, snpfname)
      call err_infplate(n, pos, cf, t, error(3))
    else if ((t-dt/2<tfinish*2/3).and.(tfinish*2/3<t+dt/2)) then
      call plot_simple(n, pos, snpfname)
      call plot_simple(n, cf, snpfname)
      call err_infplate(n, pos, cf, t, error(4))
    end if

    t = t + dt
    iter = iter + 1
  end do

  write (*, *) t - dt
  call print_output(n, t, pos, vel, acc, mas, den, h, prs, ieu, cf)
  call plot_simple(n, pos, snpfname)
  call plot_simple(n, cf, snpfname)
  call err_infplate(n, pos, cf, t, error(5))
  call plot_simple(5, error, errfname)

end program main
