program main
  use setup
  use internal
  use printer

  integer, parameter :: n = 100
  integer, parameter :: nbnd = 5
  integer, parameter :: nmax = n + 2 * nbnd

  real, parameter    :: xmin = 0.
  real, parameter    :: xmax = 1.
  real, parameter    :: init_rho = 1.
  real, parameter    :: init_smlen = 1.2
  real, parameter    :: speeedOfSound = 1.

  real :: position(nmax), velocity(nmax), density(nmax), slength(nmax)
  real :: pressure(nmax), mass(nmax), acceleration(nmax)
  real :: dt, t, dtout, ltout, kenergy
  real :: a(nmax), p(nmax), v(nmax)

  call periodic_ic(xmin, xmax, nmax, nbnd, init_rho, init_smlen, position, mass, velocity, density, slength)

  t = 0.
  dt = 0.001
  do i = 1, 10
    call output(nmax, 1. * i, position, velocity, acceleration, mass, density, slength)
    p(:) = position(:)
    v(:) = velocity(:)
    a(:) = acceleration(:)
    position(:) = p(:) + dt * v(:) + 0.5 * dt * dt * a(:)
    velocity(:) = v(:) + dt * a(:)
    call derivs(nmax, nbnd,position, mass, density, slength, pressure, acceleration, speeedOfSound)
    velocity(:) = velocity(:) + 0.5 * dt * (acceleration(:) - a(:))
    t = t + dt
  end do
  t = 0.
  dt = 0.001
  dtout = 0.5
  ltout = 0.
  ! do while (t <= 5.)
  !   if (t >= ltout) then
  !     write (*, *) t
  !     call output(nmax, t, position, velocity, acceleration, mass, density, slength)
  !     ltout = ltout + dtout
  !   end if
  !   p(:) = position(:)
  !   v(:) = velocity(:)
  !   a(:) = acceleration(:)
  !   position(:) = p(:) + dt * v(:) + 0.5 * dt * dt * a(:)
  !   velocity(:) = v(:) + dt * a(:)
  !   call derivs(nmax, nbnd,position, mass, density, slength, pressure, acceleration, speeedOfSound)
  !   velocity(:) = velocity(:) + 0.5 * dt * (acceleration(:) - a(:))
  !   call get_kinetic_energy(nmax, nbnd, mass, velocity, kenergy)
  !   call plot_simple(t, kenergy, 'energy.dat')
  !   t = t + dt
  ! end do
end program main
