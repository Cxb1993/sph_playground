program main
  use BC,               only: bcdestroy => destroy,&
                              realpartnumb
  use setup,            only: setupV2
  use state,            only: get_tasktype,&
                              ginitvar,&
                              setpartnum,&
                              getdiffisotropic, &
                              getdiffconductivity, &
                              getmhdmagneticpressure
  use iterator,         only: iterate
  use printer,          only: Output, AppendLine
  use errcalc,          only: err_diff_laplace, &
                              err_diff_graddiv, &
                              err_sinxet, &
                              err_shockTube => shockTube, &
                              err_soundwave_d => soundwaveperturbation_density, &
                              err_soundwave_v => soundwaveperturbation_velocity, &
                              err_diff_artvisc => diff_artvisc, &
                              err_alfvenwave => alfvenwave
  use args,             only: fillargs
  use errteylor,        only: etlaplace => laplace,&
                              etgraddiv => graddiv
  use circuit1,         only: c1_init, &
                              c1destroy => destroy
  use timing,           only: printTimes,&
                              tinit => init, &
                              timedestroy => destroy
  use neighboursearch,  only: getNeibNumbers,&
                              neibDestroy => destroy,&
                              getNeibListL1
  use arrayresize,      only: resize
  use timing,           only: addTime

  use preruncheck
  use const
  ! use cltest,           only: runcltest

  implicit none

  real, allocatable, dimension(:,:) :: &
    tr, tv, ta, tdb
  real, allocatable, dimension(:) :: &
    tdh, err, sqerr, tdu, tddt, &
    result
  real, allocatable, dimension(:,:) :: &
    store
  real :: &
    dt, t, dtout, ltout, tfinish, npic,&
    pspc1, gamma, sk, cv = 1., &
    mhdmuzero, difcond
    ! , chi(81)
  character (len=100) :: &
    itype, errfname, ktype, dtype
  integer :: &
    n, dim, iter, s_tt, nusedl1, nusedl2, printlen, silent,&
    ivt, stopiter, resol, difiso

  integer(8) :: &
    tprint

  ! call runcltest()
  print *, '##############################################'
  print *, '#####'

  call fillargs(dim, pspc1, resol,&
                itype, ktype, dtype, errfname, dtout, npic, tfinish, sk, silent)

  call setupV2(n, sk, gamma, cv, pspc1, resol, store)

  call setpartnum(n)
  call get_tasktype(s_tt)
  call ginitvar(ivt)
  call getdiffisotropic(difiso)
  call getdiffconductivity(difcond)
  call getmhdmagneticpressure(mhdmuzero)

  select case(s_tt)
  case (eeq_hydro, eeq_diffusion, eeq_magnetohydro, eeq_magnetohydrodiffusion)
  ! case (5, 6, 7, 8, 10)
    ! 'diff-laplass' ! 'diff-graddiv' ! diff-artvisc
    ! call set_stepping(10**dim)
  case default
    print *, 'Particle turn-off was not set main.f90: line 75.'
    stop
  end select

  print *, '#####'
  print *, '##############################################'

  allocate(result(100))
  result(:) = 0.
  result(1) = pspc1
  result(2) = realpartnumb
  n = realpartnumb

  t = 0.
  dt = 0.
  ltout = 0.
  iter = 0.
  allocate(tr(3,n))
  allocate(tv(3,n))
  allocate(ta(3,n))
  allocate(tdb(3,n))
  allocate(tdu(n))
  allocate(tdh(n))
  allocate(tddt(n))
  allocate(err(1:n))
  allocate(sqerr(1:n))

  err(:) = 0.

  call tinit()
  call c1_init(n)
  stopiter = 0
  print *, "# Finish time = ", tfinish

  if (silent == 0) then
    print *, '# Initial setup printed'
    call Output(t, store, sqerr)
  end if
  call checkVarsReady(s_tt, store)
  call iterate(n, sk, gamma, store)
  select case(s_tt)
  case(eeq_hydro, eeq_magnetohydro, eeq_diffusion, eeq_magnetohydrodiffusion)
    ! stopiter = 1
  ! case(7, 8, 9)
  !   ! who knows....
  ! case(5, 6, 10)
  !   ! 'diff-laplace' ! 'diff-graddiv' ! diff-artvisc
  !   stopiter = 1
  case default
    print *, 'EQS type was not sen in while interuption main.f90: line 139.'
    stop
  end select

  do while ((t < tfinish + eps0).and.(stopiter==0))
    ! call Output(t, ptype, pos, vel, acc, mas, den, h, prs, iu, cf, sqerr)
    select case(s_tt)
    case(eeq_hydro)
      dt = .3 * minval(store(es_h,1:n)) / maxval(store(es_c,1:n))
    case(eeq_magnetohydro)
      dt = .1 * minval(store(es_h,1:n)) / maxval(store(es_c,1:n))
    case(eeq_diffusion)
      dt = .1 * minval(store(es_den,1:n)) &
              * minval(store(es_c,1:n)) &
              * minval(store(es_h,1:n)) ** 2 &
              / merge(difcond, maxval(store(es_bx:es_bz,1:n)), difiso == 1)
    case (eeq_magnetohydrodiffusion)
      ! dt = .1 * mhd_magneticconstant * minval(den(1:realpartnumb)) * minval(c(1:realpartnumb)) * &
      !           minval(h(1:realpartnumb)) ** 2 / &
      !           merge(maxval(kcf(:,1,1:realpartnumb)), 1., maxval(kcf(:,1,1:realpartnumb))>0)
      dt = .1 * mhdmuzero &
              * minval(store(es_den,1:n)) &
              * minval(store(es_c,1:n)) &
              * minval(store(es_h,1:n)) ** 2 &
              / merge(difcond, maxval(store(es_bx:es_bz,1:n)), difiso == 1)
    ! case (5,6,7,8)
    !   ! 'diff-laplass'      ! 'diff-graddiv'
    !   stopiter = 1
    !   dt = 0.
    case default
      print *, 'Task type time increment was not set main.f90: line 150.'
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
      print *, "#", iter, "t=", t, "dt=", dt, 'min h=', minval(store(es_h,1:n)), 'max h=', maxval(store(es_h,1:n))
      if ( silent == 0) then
        call Output(t, store, err)
      end if
      ltout = ltout + dtout
    end if

    tr(:,:) = store(es_rx:es_rz,1:n)
    tv(:,:) = store(es_vx:es_vz,1:n)
    ta(:,:) = store(es_ax:es_az,1:n)
    tdu(:) = store(es_du,1:n)
    tdh(:) = store(es_dh,1:n)
    tddt(:) = store(es_ddt,1:n)
    tdb(:,:) = store(es_dbx:es_dbz,1:n)

    store(es_rx:es_rz,1:n) = store(es_rx:es_rz,1:n)  + &
        dt * tv(:,:) + 0.5 * dt * dt * ta(:,:)

    store(es_vx:es_vz,1:n) = tv(:,:)  + dt * ta(:,:)
    store(es_u,1:n) = store(es_u,1:n) + dt * tdu(:)
    store(es_h,1:n) = store(es_h,1:n) + dt * tdh(:)
    store(es_t,1:n) = store(es_t,1:n) + dt * tddt(:)
    store(es_bx:es_bz,1:n) = store(es_bx:es_bz,1:n) + dt * tdb(:,:)

    call iterate(n, sk, gamma, store)

    store(es_vx:es_vz,1:n) = &
      store(es_vx:es_vz,1:n) + 0.5*dt*(store(es_ax:es_az,1:n) - ta(:,:))
    store(es_u,1:n) = &
      store(es_u,1:n) + 0.5*dt*(store(es_du,1:n) - tdu(:))
    store(es_h,1:n) = &
      store(es_h,1:n) + 0.5*dt*(store(es_dh,1:n) - tdh(:))
    store(es_t,1:n) = &
      store(es_t,1:n) + 0.5*dt*(store(es_ddt,1:n) - tddt(:))
    store(es_bx:es_bz,1:n) = &
      store(es_bx:es_bz,1:n) + 0.5*dt*(store(es_dbx:es_dbz,1:n) - tdb(:,:))

    t = t + dt
    iter = iter + 1
  end do
  ! print*, 11111
  !----------------------------------------!
  !         l2 error calc evaluatopn       !
  !----------------------------------------!
  ! select case(s_tt)
  ! case(1, 2, 3, 4)
  !   ! mooved to ivt check below
  !   ! the rest should be mooved as well
  ! case(7, 8)
  !   ! 'hydroshock' ! chi-laplace ! 'infslb'
  !   ! there was empty
  ! case(5)
  !   ! 'diff-laplace'
  !   call err_diff_laplace(pos, acc, err)
  ! case(6)
  !   ! 'diff-graddiv'
  !   call err_diff_graddiv(ptype, pos, acc, err)
  ! case(10)
  !   ! diff-artvisc
  !   call err_diff_artvisc(pos, acc, err)
  ! case default
  !   print *, 'Task type was not set in l2 error evaluation main.f90: line 229.'
  !   stop
  ! end select

  select case(ivt)
  case (ett_pulse, ett_OTvortex)
    printlen = 5
  case (ett_sin3)
    call err_sinxet(store, t, err)
    printlen = 5
  case (ett_soundwave)
    ! call err_soundwave_v(pos, vel, t, err)
    call err_soundwave_d(store, t, err)
    printlen = 5
  case (ett_hydroshock)
    call err_shockTube(store, t, err)
    printlen = 5
  case (ett_alfvenwave)
    call err_alfvenwave(store, t, err)
    printlen = 5
  case default
    print *, ivt, 'Task type was not sen in l2 error evaluation main:267.'
    stop
  end select

  call getNeibNumbers(nusedl1, nusedl2)
  result(3) = merge(sqrt(sum(err)/nusedl1), 0., nusedl1 > 0)
  !----------------------------------------!
  !          teylor error evaluation       !
  !----------------------------------------!

  if (nusedl1 /= 0) then
    sqerr(:) = sqrt(err(:))
    print *, "#", iter, "t=", t
    if (silent == 0) then
      call Output(t, store, sqerr)
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

  call neibDestroy()
  call c1destroy()
  call timedestroy()
  call bcdestroy()
end program

subroutine set_stepping(i)
  use neighboursearch, only: stnb  => setStepsize

  integer, intent(in) :: i

  call stnb(i)
  print *, '# #      step.size:', i
end subroutine
