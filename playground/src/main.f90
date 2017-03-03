program main
  use BC
  use IC,        only: setupIC
  use kernel,    only: get_tasktype
  use iterator,  only: iterate
  use printer
  use errcalc,   only: err_init,&
                       err_diff_laplace,&
                       err_diff_graddiv,&
                       err_sinxet
  use args,      only: fillargs
  use errteylor, only: etlaplace => laplace
  use circuit1,  only: c1_init

  implicit none

  real, allocatable, dimension(:,:)  :: p, v, a
  real, allocatable, dimension(:,:)  :: pos, vel, acc
  real, allocatable, dimension(:)    :: den, prs, mas, iu, du, om, c, h, dh, &
                                        cf, dcf, kcf, tdu, tdh, tcf, err, sqerr, result
  integer, allocatable, dimension(:) :: ptype

  real                               :: dt, t, dtout, ltout, tfinish, npic,&
                                        pspc1, pspc2, gamma,&
                                        sk, chi(9)
  real                :: cv = 1.
  character (len=40)  :: itype, errfname, ktype
  integer             :: n, dim, iter, tt, nused,i

  print *, '##############################################'
  print *, '#####'
  call fillargs(dim, pspc1, pspc2,&
                itype, ktype, errfname, dtout, npic, tfinish, sk)

  call setupIC(n, sk, gamma, cv, pspc1, pspc2, pos, vel, acc, &
                mas, den, h, prs, iu, du, cf, kcf, dcf, ptype)

  call set_stepping(10**dim)
  print *, '#####'
  print *, '##############################################'

  call get_tasktype(tt)
  allocate(result(19))
  result(1) = pspc1
  result(2) = n

  t = 0.
  dt = 0.
  ltout = 0.
  iter = 0.
  nused = 0
  p = pos
  v = vel
  a = acc
  allocate(err(1:n))
  allocate(sqerr(1:n))
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

 print *, "Finish time = ", tfinish
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
      call err_sinxet(ptype, cf, t, err, nused)
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
      call print_output(t, ptype, pos, vel, acc, mas, den, h, prs, iu, cf, err)
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
    ! print *, 999
    call iterate(n, sk, gamma, &
                ptype, pos, vel, acc, &
                mas, den, h, dh, om, prs, c, iu, du, &
                cf, dcf, kcf)
    ! print *, maxval(abs(du))
    ! print *, 0, 2
    ! print *, maxval(cf)
    vel(:,:) = vel(:,:) + 0.5 * dt * (acc(:,:) - a(:,:))
    iu(:)    = iu(:)    + 0.5 * dt * (du(:) - tdu(:))
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
  end do
  ! print *, 999

  if (t <= .0) then
    select case(tt) ! l2 error calc evaluatopn
    case(1, 2, 7)
      ! 'hydroshock' ! chi-laplace ! 'infslb'
    case(3)
      ! 'hc-sinx'
      call err_sinxet(ptype, cf, t, err, nused)
    case(5)
      ! 'diff-laplace'
      call err_diff_laplace(ptype, pos, acc, err, nused)
    case(6)
      ! 'diff-graddiv'
      call err_diff_graddiv(ptype, pos, acc, err, nused)
    case default
      print *, 'Task type was not sen in error evaluation main.f90'
      stop
    end select

    select case(tt) ! teylor error evaluation
    case(7)
      ! 'hydroshock' ! chi-laplace ! 'infslb'
      call etlaplace(int(n/2), pos, mas, den, h, chi)
    case(1, 2, 3, 5, 6)
      ! 'hc-sinx' ! 'diff-laplace' ! 'diff-graddiv'
    case default
      print *, 'Task type was not sen in error evaluation main.f90'
      stop
    end select
    if (nused /= 0) then
      ! call calcDaigonal2ndErrTerms(ptype, pos, mas, den, h, chi)
      sqerr(:) = sqrt(err(:))
      call print_output(t, ptype, pos, vel, acc, mas, den, h, prs, iu, cf, sqerr)
    end if
    result(3) = merge(sqrt(sum(err)/nused), 0., nused>0)
    result(4) = 0. ! result(4) = merge(sum(chi)/nused/dim, 0., nused>0)
    result(5) = sk
    result(6:14) = chi(1:9)
    call print_appendline(14, result, errfname)
  else
    result(9) = sqrt(sum(err)/nused)
    ! call calcDaigonal2ndErrTerms(ptype, pos, mas, den, h, chi)
    ! result(10) = sum(chi)/n/dim
    result(10) = 0.
    result(11) = t
    call print_appendline(11, result, errfname)
  end if
end program main

subroutine set_stepping(i)
  use errcalc,         only: sterr => setStepsize
  use circuit2,        only: stc2  => setStepsize
  use neighboursearch, only: stnb  => setStepsize

  integer, intent(in) :: i

  call sterr(i)
  call stc2(i)
  call stnb(i)

  print *, '# #   step.size:', i
end subroutine set_stepping