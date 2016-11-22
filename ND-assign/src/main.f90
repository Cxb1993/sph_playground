program main
  use BC
  use IC
  use iterator, only:iterate
  use printer
  use err_calc
  use args,     only:fillargs
  use bias,     only:ex22
  implicit none

  real                :: sk = 1.2
  integer, parameter  :: nmax = 400000
  real                :: cv = 1.
  real                :: gamma = 1.4
  integer             :: n, dim
  character (len=40)  :: itype, errfname, ktype

  real, allocatable, dimension(:,:) :: position, velocity, acceleration
  real :: density(nmax), slength(nmax), pressure(nmax), mass(nmax), ienergy(nmax), dienergy(nmax), omega(nmax)
  real :: coupledfield(nmax), kcoupledfield(nmax), dcoupledfield(nmax)

  real                              :: dt, t, dtout, ltout, tfinish, error(11), npic, &
                                       pspc1, pspc2, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2
  integer                           :: finish, iter
  real, allocatable, dimension(:,:) :: p, v, a
  real, allocatable, dimension(:,:) :: pos, vel, acc
  real, allocatable, dimension(:)   :: den, prs, mas, ieu, diu, o, du, c, h, dh, tdh
  real, allocatable, dimension(:)   :: cf, tcf, dcf, kcf, err, ex

  allocate(position(3,nmax))
  allocate(velocity(3,nmax))
  allocate(acceleration(3,nmax))

  print *, '#######################'
  print *, '##'
  call fillargs(dim, pspc1, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc2,&
                itype, ktype, errfname, dtout, npic, tfinish)

  call setup(itype, ktype, dim, nmax, n, sk, gamma, cv, &
                pspc1, pspc2, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, &
                position, velocity, acceleration, &
                mass, density, slength, pressure, ienergy, coupledfield, kcoupledfield, dcoupledfield)
  print *, '##'
  print *, '#######################'
  read *

  error(1) = pspc1
  error(2) = n

  t = 0.
  dt = 0.
  ltout = 0.
  finish = 3000.

  pos = position(:,1:n)
  p   = pos(:,:)
  vel = velocity(:,1:n)
  v   = vel(:,:)
  acc = acceleration(:,1:n)
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
  allocate(ex(1:n))


  do while (t <= tfinish)
    select case(itype)
    case('hydroshock')
      dt = .3 * minval(h) / maxval(c)
    case('infslb', 'hc-sinx')
      dt = .144 * minval(den) * minval(c) * minval(h) ** 2 / maxval(kcf)
    case default
      print *, 'DT not set: ', itype
      stop
    end select
    ! print *, "dt: ", dt, minval(den), minval(c), minval(h), maxval(kcf)
    ! read *
    if (t >= ltout) then
      print *, iter, t
      select case(itype)
      case('hydroshock')
        call err_sinxet(n, pos, cf, t, err)
      case('infslb')
        call err_sinxet(n, pos, cf, t, err)
      case('hc-sinx')
        call err_sinxet(n, pos, cf, t, err)
      case default
        print *, 'Main: type not set ', itype
        stop
      end select
      call print_output(n, t, pos, vel, acc, mas, den, h, prs, ieu, cf, sqrt(err))
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

    call iterate(n, sk, gamma, &
                pos, vel, acc, &
                mas, den, h, dh, o, prs, c, ieu, diu, &
                cf, dcf, kcf)

    vel(:,:) = vel(:,:) + 0.5 * dt * (acc(:,:) - a(:,:))
    ieu(:)   = ieu(:)   + 0.5 * dt * (diu(:) - du(:))
    h(:)     = h(:)     + 0.5 * dt * (dh(:) - tdh(:))

    select case(itype)
    case('hydroshock')
      cf(:) = cf(:)     + 0.5 * dt * (dcf(:) - tcf(:))
    case('infslb', 'hc-sinx')
      cf(:) = ieu(:) / cv
    case default
      print *, 'Coupled field integration not set: ', itype
      stop
    end select

    t = t + dt
    iter = iter + 1

    if ((t-dt/2<tfinish*1/3).and.(tfinish*1/3<t+dt/2)) then
      error(3) = sqrt(sum(err)/n)
      call ex22(n, mas, den, pos, h, ex)
      call periodic1(ex, 1)
      error(4) = sum(ex)/n
      error(5) = t
    else if ((t-dt/2<tfinish*2/3).and.(tfinish*2/3<t+dt/2)) then
      error(6) = sqrt(sum(err)/n)
      call ex22(n, mas, den, pos, h, ex)
      error(7) = sum(ex)/n
      error(8) = t
    end if
  end do

  error(9) = sqrt(sum(err)/n)
  call ex22(n, mas, den, pos, h, ex)
  error(10) = sum(ex)/n
  error(11) = t
  call plot_simple(11, error, errfname)
end program main
