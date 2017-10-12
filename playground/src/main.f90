program main
  use BC,               only: bcdestroy => destroy
  use IC,               only: setupIC
  use state,            only: get_tasktype,&
                              ginitvar
  use iterator,         only: iterate
  use printer,          only: Output, AppendLine
  use errcalc,          only: err_diff_laplace, &
                              err_diff_graddiv, &
                              err_sinxet, &
                              err_shockTube => shockTube, &
                              err_soundwave_d => soundwaveperturbation_density, &
                              err_soundwave_v => soundwaveperturbation_velocity
  use args,             only: fillargs
  use errteylor,        only: etlaplace => laplace,&
                              etgraddiv => graddiv
  use circuit1,         only: c1_init, &
                              c1destroy => destroy
  use timing,           only: printTimes,&
                              tinit => init, &
                              timedestroy => destroy
  use neighboursearch,  only: getNeibNumbers,&
                              neibDestroy => destroy
  use arrayresize,      only: resize
  use map,              only: appendmap,&
                              mapdestroy => destroy
  use const

  use cltest,           only: runcltest

  implicit none

  real, allocatable, dimension(:,:,:) :: dfdx, kcf
  real, allocatable, dimension(:,:)   :: p, v, a, pos, vel, acc, cf, dcf, tcf
  real, allocatable, dimension(:)     :: den, prs, mas, iu, du, om, c, h, dh, &
                                         tdh, err, sqerr, tdu,&
                                         result
  integer, allocatable, dimension(:)  :: ptype

  real                                :: dt, t, dtout, ltout, tfinish, npic,&
                                         pspc1, pspc2, gamma,&
                                         sk, chi(81), cv = 1.

  character (len=100)  :: itype, errfname, ktype, dtype
  integer             :: n, dim, iter, tt, nusedl1, nusedl2, printlen, silent, ivt, stopiter

  integer(8)          :: tprint

  ! call runcltest()

  print *, '##############################################'
  print *, '#####'

  call fillargs(dim, pspc1, pspc2,&
                itype, ktype, dtype, errfname, dtout, npic, tfinish, sk, silent)

  call setupIC(n, sk, gamma, cv, pspc1, pspc2, pos, vel, acc, &
                mas, den, h, prs, iu, du, cf, kcf, dcf, ptype)
  call get_tasktype(tt)
  call ginitvar(ivt)

  select case(tt)
  case (1, 2, 3, 4, 9)
  case (5, 6, 7, 8)
    ! 'diff-laplass'      ! 'diff-graddiv'
    call set_stepping(10**dim)
  case default
    print *, 'Particle turn-off was not set main.f90: line 56.'
    stop
  end select

  print *, '#####'
  print *, '##############################################'

  allocate(result(100))
  result(:) = 0.
  result(1) = pspc1
  result(2) = n

  t = 0.
  dt = 0.
  ltout = 0.
  iter = 0.
  allocate(p(3,n))
  allocate(v(3,n))
  allocate(a(3,n))
  p = pos
  v = vel
  a = acc
  allocate(err(1:n))
  err = 0.
  allocate(sqerr(1:n))
  allocate(c(n))
  c(:) = cv
  allocate(tdu(n))
  allocate(dh(n))
  dh(:) = 0
  allocate(tdh(n))
  allocate(tcf(3,n))
  allocate(om(n))
  if ((dtype == 'symm') .or. (ktype == '2nw')) then
    allocate(dfdx(3,3,n))
  end if

  call tinit()
  call c1_init(n)
  stopiter = 0
  print *, "Finish time = ", tfinish

  call iterate(n, sk, gamma, &
              ptype, pos, vel, acc, &
              mas, den, h, dh, om, prs, c, iu, du, &
              cf, dcf, kcf, dfdx)

  select case(tt)
  case(1, 2, 3, 7, 8, 9)
    ! 'hydroshock' ! 'heatconduction' ! chi-laplace ! 'infslb'
  case(5, 6)
    ! 'diff-laplace' ! 'diff-graddiv'
    stopiter = 1
  case default
    print *, 'Task type was not sen in l2 error evaluation main.f90: line 182'
    stop
  end select

  do while ((t < tfinish + epsilon(0.)).and.(stopiter==0))
    select case(tt)
    case(1, 9)
      ! 'hydroshock'
      dt = .3 * minval(h) / maxval(c)
    case(2)
      ! 'infslb'
      dt = .1 * minval(den) * minval(c) * minval(h) ** 2 / maxval(kcf)
    case(3)
      ! 'hc-sinx'
      dt = .1 * minval(den) * minval(c) * minval(h) ** 2 / maxval(kcf)
    case(4)
      ! 'photoevaporation' 'pheva'
      dt = .3e-3 * minval(h)**2 / maxval(c)**2 / maxval(kcf) / maxval(cf)
    case (5,6,7,8)
      ! 'diff-laplass'      ! 'diff-graddiv'
      stopiter = 1
      dt = 0.
    case default
      print *, 'Task type time increment was not set main.f90: line 115.'
      stop
    end select
    if (t + dt > tfinish) then
      dt = tfinish - t
      stopiter = 1
    end if
    ! print *, 0, 0
    if (t >= ltout) then
      ! print*, maxval(kcf), minval(den), minval(c), minval(h)
      ! read*
      print *, iter, t, dt, sum(iu)
      if ( silent == 0) then
        call Output(t, ptype, pos, vel, acc, mas, den, h, prs, iu, cf, err)
      end if
      ltout = ltout + dtout
    end if

    p(:,:) = pos(:,:)
    v(:,:) = vel(:,:)
    a(:,:) = acc(:,:)
    tdu(:) = du(:)
    tdh(:) = dh(:)
    tcf(:,:) = dcf(:,:)
    pos(:,:) = p(:,:)  + dt * v(:,:) + 0.5 * dt * dt * a(:,:)
    vel(:,:) = v(:,:)  + dt * a(:,:)
    iu(:)    = iu(:)   + dt * du(:)
    h(:)     = h(:)    + dt * dh(:)
    cf(:,:)  = cf(:,:) + dt * dcf(:,:)
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
    vel(:,:) = vel(:,:) + 0.5 * dt * (acc(:,:)  - a(:,:))
    iu(:)    = iu(:)    + 0.5 * dt * (du(:)     - tdu(:))
    h(:)     = h(:)     + 0.5 * dt * (dh(:)     - tdh(:))
    cf(:,:)  = cf(:,:)  + 0.5 * dt * (dcf(:,:)  - tcf(:,:))

    t = t + dt
    iter = iter + 1
  end do
  ! print*, 11111
  !----------------------------------------!
  !         l2 error calc evaluatopn       !
  !----------------------------------------!
  select case(tt)
  case(1)
    call err_shockTube(ptype, pos, den, t, err)
  case(2, 7, 8)
    ! 'hydroshock' ! chi-laplace ! 'infslb'
  case(3)
    ! 'heatconduction'
    select case(ivt)
    case(1)
      call err_sinxet(pos, cf, t, err)
    end select
  case(5)
    ! 'diff-laplace'
    call err_diff_laplace(pos, acc, err)
  case(6)
    ! 'diff-graddiv'
    call err_diff_graddiv(ptype, pos, acc, err)
  case(9)
    ! call err_soundwave_v(pos, vel, t, err)
    call err_soundwave_d(pos, den, t, err)
    ! print*, err
  case default
    print *, 'Task type was not sen in l2 error evaluation main.f90: line 182'
    stop
  end select
  call getNeibNumbers(nusedl1, nusedl2)
  result(3) = merge(sqrt(sum(err)/nusedl1), 0., nusedl1 > 0)
  !----------------------------------------!
  !          teylor error evaluation       !
  !----------------------------------------!
  select case(tt)
  case(5, 7)
    ! diff-laplace ! chi-laplace
    call etlaplace(pos, mas, den, h, chi)
    result(4) = sum(chi(1:9))/dim
    result(6:14) = chi(1:9)
    printlen = 14
    ! print(result(4))
  case(6, 8)
    ! diff-graddiv ! chi-graddiv
    call etgraddiv(pos, mas, den, h, chi)
    result(4) = sum(chi)/dim/(2*dim - 1)
    result(6:86) = chi(1:81)
    printlen = 86
  case(1, 2, 3, 9)
    ! 'hc-sinx' ! 'diff-graddiv'
    printlen = 5
  case default
    print *, 'Task type was not sen in taylor error evaluation main.f90: line 208'
    stop
  end select
  if (nusedl1 /= 0) then
    sqerr(:) = sqrt(err(:))
    print *, iter, t, 0.
    if (silent == 0) then
      call Output(t, ptype, pos, vel, acc, mas, den, h, prs, iu, cf, sqerr)
    end if
  end if
  result(5) = sk
  call resize(result, printlen, printlen)
  call AppendLine(result, errfname, tprint)

  call printTimes()
  print *, '#####  Results:'
  write(*, "(A, F20.10)") " # #   l2-error: ", result(3)
  write(*, "(A, F20.10)") " # #  chi-error: ", result(4)
  print *, '##############################################'

  deallocate(result)
  deallocate(ptype)
  deallocate(pos)
  deallocate(vel)
  deallocate(acc)
  deallocate(mas)
  deallocate(den)
  deallocate(h)
  deallocate(prs)
  deallocate(iu)
  deallocate(du)
  deallocate(cf)
  deallocate(kcf)
  deallocate(dcf)
  deallocate(err)
  deallocate(sqerr)
  deallocate(c)
  deallocate(tdu)
  deallocate(dh)
  deallocate(tdh)
  deallocate(tcf)
  deallocate(om)
  deallocate(p)
  deallocate(v)
  deallocate(a)
  if ( ktype == '2nw' ) then
    deallocate(dfdx)
  end if
  call neibDestroy()
  call c1destroy()
  call timedestroy()
  call bcdestroy()
  call mapdestroy()
end program

subroutine set_stepping(i)
  use neighboursearch, only: stnb  => setStepsize

  integer, intent(in) :: i

  call stnb(i)
  print *, '# #      step.size:', i
end subroutine
