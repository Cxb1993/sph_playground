program main
  use BC
  use IC_hydro_shock
  use IC_heat_cond
  use internal
  use printer
  use kernel

  implicit none

  character (len=40) :: ttype
  character (len=1)  :: arg

  real                :: sk = 1.2
  integer, parameter  :: nmax = 400000
  real                :: speedOfSound = 1.
  real                :: gamma = 1.4
  integer             :: n, dim

  real, allocatable, dimension(:,:) :: position, velocity, acceleration
  ! real :: position(nmax,3), velocity(nmax,3), acceleration(nmax,3)
  real :: density(nmax), slength(nmax), pressure(nmax), mass(nmax), ienergy(nmax), dienergy(nmax), omega(nmax)
  real :: coupledfield(nmax), kcoupledfield(nmax)

  real :: dt, t, dtout, ltout, tfinish
  real, allocatable, dimension(:,:) :: p, v, a
  real, allocatable, dimension(:,:) :: pos, vel, acc
  real, allocatable, dimension(:)   :: den, prs, mas, ieu, diu, o, du, c, h, dh, tdh
  real, allocatable, dimension(:)   :: cf, tcf, dcf, kcf

  allocate(position(nmax,3))
  allocate(velocity(nmax,3))
  allocate(acceleration(nmax,3))

  call get_command_argument(1, arg)
  read(arg(:), fmt="(i5)") dim
  print *, "#       dim:", dim

  call get_command_argument(2, ttype)
  print *, "# task type:   ", ttype

  ! call periodic_ic(xmin, xmax, nmax, nbnd, init_rho, sk, &
  !                  position, mass, velocity, acceleration, density, slength, pressure)
  ! ptype='periodic'
  select case(ttype)
  case('hydroshock')
    call set_tasktype(ttype)
    call setup_hydro_shock(dim, nmax, n, sk, gamma, &
                  position, velocity, acceleration, &
                  mass, density, slength, pressure, ienergy)
    tfinish = 0.3
    dtout = 0.001
  case('temperhomog01')
    call set_tasktype(ttype)
    call setup_heat_cond(dim, nmax, n, sk, gamma, &
                  position, velocity, acceleration, &
                  mass, density, slength, pressure, ienergy, coupledfield, kcoupledfield)
    tfinish = 650
    dtout = 1
  end select

  pos = position(1:n,:)
  p   = pos(:,:)
  vel = velocity(1:n,:)
  v   = vel(:,:)
  acc = acceleration(1:n,:)
  a   = acc(:,:)
  den = density(1:n)
  h   = slength(1:n)
  prs = pressure(1:n)
  mas = mass(1:n)
  ieu = ienergy(1:n)
  diu = dienergy(1:n)
  du  = diu(:)
  o   = omega(1:n)
  allocate(c(1:n))
  c(:)= speedOfSound
  allocate(dh(1:n))
  dh  = 0
  tdh = dh
  cf  = coupledfield(1:n)
  tcf = coupledfield(1:n)
  allocate(dcf(1:n))
  dcf = 0.
  kcf  = kcoupledfield(1:n)

  deallocate(position)
  deallocate(velocity)
  deallocate(acceleration)

  t = 0.
  ltout = 0.
  call output(n, 0., pos, vel, acc, mas, den, h, prs, ieu, cf)
  call derivs(n, speedOfSound, sk, gamma, &
              pos, vel, acc, &
              mas, den, h, dh, o, prs, c, ieu, diu, &
              cf, dcf, kcf)
  print *, ''
  do while (t <= tfinish)
    select case(ttype)
    case('hydroshock')
      dt = .3 * minval(h) / maxval(c)
    case('temperhomog01')
      dt = .144 * minval(den) * minval(c) * minval(h) ** 2 / maxval(kcf)
    end select
    ! print *, "dt: ", dt, minval(den), minval(c), minval(h), maxval(kcf)
    ! read *
    if (t >= ltout) then
      write (*, *) t
      call output(n, t, pos, vel, acc, mas, den, h, prs, ieu, cf)
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

    call derivs(n, speedOfSound, sk, gamma, &
                pos, vel, acc, &
                mas, den, h, dh, o, prs, c, ieu, diu, &
                cf, dcf, kcf)

    vel(:,:) = vel(:,:) + 0.5 * dt * (acc(:,:) - a(:,:))
    ieu(:)   = ieu(:)   + 0.5 * dt * (diu(:) - du(:))
    h(:)     = h(:)     + 0.5 * dt * (dh(:) - tdh(:))
    ! cf(:)    = cf(:)    + 0.5 * dt * (dcf(:) - tcf(:))
    ! cf(:)    = ieu(:) / c(:)

    t = t + dt
  end do
  write (*, *) t - dt
  call output(n, t, pos, vel, acc, mas, den, h, prs, ieu, cf)
end program main
