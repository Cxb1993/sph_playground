program main
  use setup
  use internal
  use printer
  use kernel

  implicit none

  character (len=40) :: ptype
  character (len=1) :: arg

  real                :: sk = 1.2
  integer, parameter  :: nmax = 400000
  real                :: speedOfSound = 1.
  real                :: gamma = 1.4
  integer             :: n, dim

  real :: position(nmax,3), velocity(nmax,3), acceleration(nmax,3), density(nmax), slength(nmax)
  real :: pressure(nmax), mass(nmax), ienergy(nmax), dienergy(nmax), omega(nmax), actacc(nmax)
  real :: dustfraction(nmax)

  real :: dt, t, dtout, ltout, tfinish
  real, allocatable, dimension(:,:) :: p, v, a
  real, allocatable, dimension(:,:) :: pos, vel, acc
  real, allocatable, dimension(:)   :: den, prs, mas, ieu, diu, o, du, c, h, dh, tdh, f, eps

  call get_command_argument(1, arg)
  read(arg(:), fmt="(i5)") dim
  print *, "#      dim:", dim

  ! call periodic_ic(xmin, xmax, nmax, nbnd, init_rho, sk, &
  !                  position, mass, velocity, acceleration, density, slength, pressure)
  ! ptype='periodic'
  call shock_ic(dim, nmax, n, sk, gamma, &
                position, velocity, acceleration, &
                mass, density, slength, pressure, ienergy, actacc, dustfraction)
  ptype='shock_fixed'
  
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
  tdh = dh
  f   = actacc(1:n)
  eps = dustfraction(1:n)

  tfinish = 0.3
  t = 0.
  dtout = 0.001
  ltout = 0.
  call output(n, 0., pos, vel, acc, mas, den, h, prs, ieu)
  call derivs(ptype, n, speedOfSound, sk, gamma, &
              pos, vel, acc, &
              mas, den, h, dh, o, prs, c, ieu, diu, f, eps)

  print *, ''
  do while (t <= tfinish)
    dt = 0.3 * minval(h) / maxval(c)
    if (t >= ltout) then
      write (*, *) t
      call output(n, t, pos, vel, acc, mas, den, h, prs, ieu)
      ltout = ltout + dtout
    end if
    p(:,:) = pos(:,:)
    v(:,:) = vel(:,:)
    a(:,:) = acc(:,:)
    du(:)  = diu(:)
    tdh(:) = dh(:)

    pos(:,:) = p(:,:) + dt * v(:,:) + 0.5 * dt * dt * a(:,:)
    vel(:,:) = v(:,:) + dt * a(:,:)
    ieu(:)   = ieu(:) + dt * diu(:)
    h(:)     = h(:)   + dt * dh(:)

    call derivs(ptype, n, speedOfSound, sk, gamma, &
                pos, vel, acc, &
                mas, den, h, dh, o, prs, c, ieu, diu, f, eps)

    vel(:,:) = vel(:,:) + 0.5 * dt * (acc(:,:) - a(:,:))
    ieu(:)   = ieu(:)   + 0.5 * dt * (diu(:) - du(:))
    h(:)     = h(:)     + 0.5 * dt * (dh(:) - tdh(:))

    t = t + dt
  end do
  write (*, *) t - dt
  call output(n, t, pos, vel, acc, mas, den, h, prs, ieu)
end program main
