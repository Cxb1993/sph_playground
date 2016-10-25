program main
  use BC
  use IC
  use internal
  use printer
  use kernel
  use err_calc

  implicit none

  character (len=40) :: itype
  character (len=1)  :: arg

  real                :: sk = 1.2
  integer, parameter  :: nmax = 400000
  real                :: cv = 1.
  real                :: gamma = 1.4
  integer             :: n, dim

  real, allocatable, dimension(:,:) :: position, velocity, acceleration
  ! real :: position(nmax,3), velocity(nmax,3), acceleration(nmax,3)
  real :: density(nmax), slength(nmax), pressure(nmax), mass(nmax), ienergy(nmax), dienergy(nmax), omega(nmax)
  real :: coupledfield(nmax), kcoupledfield(nmax), dcoupledfield(nmax)

  real :: dt, t, dtout, ltout, tfinish
  real, allocatable, dimension(:,:) :: p, v, a
  real, allocatable, dimension(:,:) :: pos, vel, acc
  real, allocatable, dimension(:)   :: den, prs, mas, ieu, diu, o, du, c, h, dh, tdh
  real, allocatable, dimension(:)   :: cf, tcf, dcf, kcf, err

  allocate(position(nmax,3))
  allocate(velocity(nmax,3))
  allocate(acceleration(nmax,3))

  call get_command_argument(1, arg)
  read(arg(:), fmt="(i5)") dim
  print *, "#       dim:", dim

  call get_command_argument(2, itype)
  print *, "# task type:   ", itype

  call setup(itype, dim, nmax, n, sk, gamma, cv, &
                position, velocity, acceleration, &
                mass, density, slength, pressure, ienergy, coupledfield, kcoupledfield, dcoupledfield)

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
  allocate(err(1:n))

  t = 0.
  ltout = 0.
  call print_output(n, 0., pos, vel, acc, mas, den, h, prs, ieu, cf, err)
  call derivs(n, sk, gamma, &
              pos, vel, acc, &
              mas, den, h, dh, o, prs, c, ieu, diu, &
              cf, dcf, kcf)
  print *, ''

  tfinish = 0.
  dt = 0.
  dtout = 0.
  select case(itype)
  case('hydroshock')
    tfinish = 0.4
  case('heatslab')
    tfinish = 7.
  end select
  do while (t <= tfinish)
    select case(itype)
    case('hydroshock')
      dt = .3 * minval(h) / maxval(c)
      dtout = 10. * dt
    case('heatslab')
      dt = .144 * minval(den) * minval(c) * minval(h) ** 2 / maxval(kcf)
      dtout = 20. * dt
    end select
    ! print *, "dt: ", dt, minval(den), minval(c), minval(h), maxval(kcf)
    ! read *
    if (t >= ltout) then
      write (*, *) t
      call print_output(n, t, pos, vel, acc, mas, den, h, prs, ieu, cf, err)
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
      call err_T0sxsyet(n, pos, cf, t, err)
    else
      cf(:) = cf(:)     + 0.5 * dt * (dcf(:) - tcf(:))
    end if

    t = t + dt
  end do
  write (*, *) t - dt
  call print_output(n, t, pos, vel, acc, mas, den, h, prs, ieu, cf, err)
end program main
