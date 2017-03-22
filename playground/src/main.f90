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
  use errteylor, only: etlaplace => laplace,&
                       etgraddiv => graddiv
  use circuit1,  only: c1_init

  implicit none

  real, allocatable, dimension(:,:,:) :: dfdx
  real, allocatable, dimension(:,:)   :: p, v, a, pos, vel, acc
  real, allocatable, dimension(:)     :: den, prs, mas, iu, du, om, c, h, dh, &
                                         cf, dcf, kcf, tdu, tdh, tcf, err, sqerr,&
                                         result
  integer, allocatable, dimension(:)  :: ptype

  real                                :: dt, t, dtout, ltout, tfinish, npic,&
                                         pspc1, pspc2, gamma,&
                                         sk, chi(81)
  real                :: cv = 1.
  character (len=40)  :: itype, errfname, ktype, dtype
  integer             :: n, dim, iter, tt, nused, printlen!, i

  print *, '##############################################'
  print *, '#####'
  call fillargs(dim, pspc1, pspc2,&
                itype, ktype, dtype, errfname, dtout, npic, tfinish, sk)

  call setupIC(n, sk, gamma, cv, pspc1, pspc2, pos, vel, acc, &
                mas, den, h, prs, iu, du, cf, kcf, dcf, ptype)

  call set_stepping(10**dim)
  ! call set_stepping(2)
  print *, '#####'
  print *, '##############################################'

  call get_tasktype(tt)
  allocate(result(100))
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
  allocate(dfdx(3,3,n))

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
                cf, dcf, kcf, dfdx)
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
  ! print*, 11111
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
    result(3) = merge(sqrt(sum(err)/nused), 0., nused>0)
    !----------------------------------------!
    !          teylor error evaluation       !
    !----------------------------------------!
    select case(tt)
    case(5)
      ! 'diff-laplace'
      call etlaplace(pos, mas, den, h, chi)
      result(4) = merge(sum(chi(1:9))/dim, 0., nused>0)
      result(6:14) = chi(1:9)
      printlen = 14
    case(6)
      ! 'diff-graddiv'
      call etgraddiv(pos, mas, den, h, chi)
      if ( dim == 1) then
        result(4) = sum(chi)
      elseif ( dim == 2 ) then
        result(4) = merge(sum(chi)/dim/3., 0., nused>0)
      elseif ( dim == 3 ) then
        result(4) = merge(sum(chi)/dim/5., 0., nused>0)
      end if
      result(6:86) = chi(1:81)
      printlen = 86
    case(7)
      ! 'hydroshock' ! chi-laplace ! 'infslb'
      call etlaplace(pos, mas, den, h, chi)
      result(4) = merge(sum(chi)/dim, 0., nused>0)
      result(6:14) = chi(1:9)
      printlen = 14
    case(1, 2, 3)
      ! 'hc-sinx' ! 'diff-graddiv'
    case default
      print *, 'Task type was not sen in error evaluation main.f90'
      stop
    end select
    if (nused /= 0) then
      sqerr(:) = sqrt(err(:))
      call print_output(t, ptype, pos, vel, acc, mas, den, h, prs, iu, cf, sqerr)
    end if
    result(5) = sk
    call print_appendline(printlen, result, errfname)
  else
    result(9) = sqrt(sum(err)/nused)
    result(10) = 0.
    result(11) = t
    call print_appendline(11, result, errfname)
  end if
  call getTime()
  print *, '##### Results:'
  write(*, "(A, F10.5)") " # #  l2-error: ", result(3)
  write(*, "(A, F10.5)") " # # chi-error: ", result(4)
  ! write(*, "(A, F10.5)") " # #     neibs: ", elapsed
  print *, '##############################################'
end program main

subroutine set_stepping(i)
  use kernel,          only: get_kerntype
  use errcalc,         only: sterr => setStepsize
  use circuit1,        only: stc1  => setStepsize
  use circuit2,        only: stc2  => setStepsize
  use neighboursearch, only: stnb  => setStepsize

  integer, intent(in) :: i
  integer :: ktp

  call get_kerntype(ktp)
  call sterr(i)
  call stc1(i)
  call stc2(i)
  call stnb(i)
  print *, '# #   step.size:', i
end subroutine set_stepping

subroutine getTime()
  use circuit1,        only: c1time => getTime
  use circuit2,        only: c2time => getTime
  use neighboursearch, only: nbtime => getTime

  real :: elapsed

  print *, '##############################################'
  print *, '#####    Time:'
  call c1time(elapsed)
  write(*, "(A, F10.5)") " # #        c1: ", elapsed
  call c2time(elapsed)
  write(*, "(A, F10.5)") " # #        c2: ", elapsed
  call nbtime(elapsed)
  write(*, "(A, F10.5)") " # #     neibs: ", elapsed
  print *, '##############################################'

end subroutine getTime
