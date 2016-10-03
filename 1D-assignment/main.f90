program main
  use setup
  use internal
  use printer
  use kernel

  character (len=40) :: ptype

  ! integer, parameter :: n = 100
  ! integer, parameter :: nbnd = 5
  ! integer, parameter :: nmax = n + 2 * nbnd
  ! real, parameter    :: xmin = 0.
  ! real, parameter    :: xmax = 1.
  ! real, parameter    :: init_rho = 1.
  ! real, parameter    :: sk = 1.2
  ! real, parameter    :: speedOfSound = 1.

  real, parameter    :: shock_xa = 0.5
  real, parameter    :: shock_xb = 0.48
  integer, parameter :: shock_na = 500
  integer, parameter :: shock_nb = 60
  integer, parameter :: shock_nl = 6
  real, parameter    :: shock_pa = 1.
  real, parameter    :: shock_pb = 0.1
  integer, parameter :: nbnd = 6
  real, parameter    :: shock_rhoa = 1.
  real, parameter    :: shock_rhob = 0.125
  real, parameter    :: sk = 1.2
  integer, parameter :: nmax = shock_na+shock_nb+1
  real, parameter    :: speedOfSound = 1.
  real, parameter    :: gamma = 1.4


  real :: position(nmax), velocity(nmax), density(nmax), ddensity(nmax), slength(nmax), omega(nmax)
  real :: pressure(nmax), c(nmax), mass(nmax), acceleration(nmax), ienergy(nmax), dienergy(nmax)
  real :: dt, t, dtout, ltout
  real :: a(nmax), p(nmax), v(nmax), du(nmax), drho(nmax)

  ! call periodic_ic(xmin, xmax, nmax, nbnd, init_rho, sk, &
  !                  position, mass, velocity, acceleration, density, slength, pressure)
  ! ptype='periodic'
  call shock_ic(shock_xa, shock_xb, shock_na, shock_nb, shock_rhoa, shock_rhob, shock_pa, shock_pb, sk, &
                position, mass, velocity, acceleration, density, slength, pressure, ienergy, gamma)
  ptype='fixed'

  tfinish = 0.3
  t = 0.
  dtout = 0.001
  ltout = 0.

  c(:) = speedOfSound

  call derivs(ptype, nmax, nbnd, &
              position, velocity, mass, density, ddensity, omega, slength, pressure, c, acceleration, &
              ienergy, dienergy, speedOfSound, sk, gamma)
  do while (t <= tfinish)
    dt = 0.3 * minval(slength) / maxval(c)
    if (t >= ltout) then
      write (*, *) t
      call output(nmax, t, position, velocity, acceleration, mass, density, slength, pressure, ienergy)
      ltout = ltout + dtout
    end if

    p(:) = position(:)
    v(:) = velocity(:)
    a(:) = acceleration(:)
    du(:) = dienergy(:)
    drho(:) = ddensity(:)
    position(:) = p(:) + dt * v(:) + 0.5 * dt * dt * a(:)
    velocity(:) = v(:) + dt * a(:)
    ienergy(:) = ienergy(:) + dt * dienergy(:)
    density(:) = density(:) + dt * ddensity(:)
    call derivs(ptype, nmax, nbnd, &
                position, velocity, mass, density, ddensity, omega, slength, pressure, c, acceleration, &
                ienergy, dienergy, speedOfSound, sk, gamma)
    velocity(:) = velocity(:) + 0.5 * dt * (acceleration(:) - a(:))
    ienergy(:) = ienergy(:) + 0.5 * dt * (dienergy(:) - du(:))
    density(:) = density(:) + 0.5 * dt * (ddensity(:) - drho(:))
  !   call get_kinetic_energy(nmax, nbnd, mass, velocity, kenergy)
  !   call plot_simple(t, kenergy, 'energy.dat')
    t = t + dt
  end do
  write (*, *) t
  call output(nmax, t, position, velocity, acceleration, mass, density, slength, pressure, ienergy)
end program main
