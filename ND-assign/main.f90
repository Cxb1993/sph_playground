program main
  use setup
  use internal
  use printer
  use kernel

  character (len=40) :: ptype
  character (len=1) :: arg
  ! integer, parameter :: n = 100
  ! integer, parameter :: nbnd = 5
  ! integer, parameter :: nmax = n + 2 * nbnd
  ! real, parameter    :: xmin = 0.
  ! real, parameter    :: xmax = 1.
  ! real, parameter    :: init_rho = 1.
  ! real, parameter    :: sk = 1.2
  ! real, parameter    :: speedOfSound = 1.

  integer             :: nb = 6
  real                :: sk = 1.2
  integer, parameter  :: nmax = 100000
  real                :: speedOfSound = 1.
  real                :: gamma = 1.4
  integer             :: n, i, dim

  real :: position(nmax,3), velocity(nmax,3), acceleration(nmax,3), density(nmax), slength(nmax)
  real :: pressure(nmax), mass(nmax), ienergy(nmax), dienergy(nmax), omega(nmax), c(nmax)
  real :: dt, t, dtout, ltout
  real :: p(nmax,3), v(nmax,3), a(nmax,3), du(nmax)

  call get_command_argument(1, arg)
  read(arg(:), fmt="(i5)") dim
  print *, "Number of dimmentions:", dim

  ! call periodic_ic(xmin, xmax, nmax, nbnd, init_rho, sk, &
  !                  position, mass, velocity, acceleration, density, slength, pressure)
  ! ptype='periodic'
  call shock_ic(dim, nmax, n, sk, gamma, &
                position, velocity, acceleration, mass, density, slength, pressure, ienergy)
  ptype='shock_fixed'
  c(:) = speedOfSound

  tfinish = 0.3
  t = 0.
  dtout = 0.001
  ltout = 0.

  call output(n, 0., position(1:n,:), velocity(1:n,:), acceleration(1:n,:), mass, density, slength, pressure, ienergy)
  call derivs(ptype, n, nb, &
              position(1:n,:), velocity(1:n,:), acceleration(1:n,:), mass, density, slength, omega, pressure, c, ienergy, &
              dienergy, speedOfSound, sk, gamma)

  do while (t <= tfinish)
    dt = 0.3 * minval(slength(1:n)) / maxval(c(1:n))
    if (t >= ltout) then
      write (*, *) t
      call output(n, t, position(1:n,:), velocity(1:n,:), acceleration(1:n,:), mass, density, slength, pressure, ienergy)
      ltout = ltout + dtout
    end if
    p(:,:) = position(:,:)
    v(:,:) = velocity(:,:)
    a(:,:) = acceleration(:,:)
    du(:) = dienergy(:)
    position(:,:) = p(:,:) + dt * v(:,:) + 0.5 * dt * dt * a(:,:)
    velocity(:,:) = v(:,:) + dt * a(:,:)
    ienergy(:) = ienergy(:) + dt * dienergy(:)
    call derivs(ptype, n, nb, &
                position(1:n,:), velocity(1:n,:), acceleration(1:n,:), mass, density, slength, omega, pressure, c, ienergy, &
                dienergy, speedOfSound, sk, gamma)
    velocity(:,:) = velocity(:,:) + 0.5 * dt * (acceleration(:,:) - a(:,:))
    ienergy(:) = ienergy(:) + 0.5 * dt * (dienergy(:) - du(:))
    t = t + dt
  end do
  write (*, *) t - dt
  call output(n, t, position(1:n,:), velocity(1:n,:), acceleration(1:n,:), mass, density, slength, pressure, ienergy)
end program main
