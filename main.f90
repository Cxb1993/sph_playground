program main
  use internal
  use setup
  use printer

  integer, parameter :: n = 100
  integer, parameter :: nbnd = 5
  integer, parameter :: nmax = n + 2 * nbnd
  real, parameter    :: xmin = 0.0
  real, parameter    :: xmax = 1.0

  real, parameter    :: init_rho = 1.0
  real, parameter    :: init_smlen = 1.2
  real, parameter    :: speeedOfSound = 1.0

  real :: position(nmax), velocity(nmax), density(nmax), slength(nmax)
  real :: pressure(nmax), mass(nmax), acceleration(nmax)
  real :: dt, t, dtout, ltout, kenergy
  real :: a(nmax), p(nmax), v(nmax)

  call periodic(xmin, xmax, nmax, nbnd, init_rho, init_smlen, position, mass, velocity, density, slength)
  call derivs(position, mass, density, slength, pressure, acceleration, nmax-1, nbnd, speeedOfSound)
  t = 0.0
  call output(position, velocity, acceleration, mass, density, slength, t, nmax-1)
  dt = 0.0001
  dtout = 1.0
  ltout = 0.0
  do while (t < 5.0)
    t = t + dt
    if (t > ltout) then
      ltout = ltout + dtout
      call output(position, velocity, acceleration, mass, density, slength, t, nmax-1)
      write (*, *) t
    end if
    p(:) = position(:)
    v(:) = velocity(:)
    a(:) = acceleration(:)
    position(:) = p(:) + dt * v(:) + 0.5 * dt * dt * a(:)
    velocity(:) = v(:) + dt * a(:)
    call derivs(position, mass, density, slength, pressure, acceleration, nmax-1, nbnd, speeedOfSound)
    velocity(:) = velocity(:) + 0.5 * dt * (acceleration(:) - a(:))
    call get_kinetic_energy(mass, velocity, kenergy, nmax-1, nbnd)
    call plot_simple(t, kenergy, 'energy.dat')
  end do
end program main
