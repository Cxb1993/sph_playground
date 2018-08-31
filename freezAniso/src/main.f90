program main
  ! use cudafor
  use BC,               only: initBorders
  use setup,            only: setupV2
  use state,            only: get_equations,&
                              ginitvar,&
                              getdiffisotropic, &
                              getdiffconductivity, &
                              getmhdmagneticpressure,&
                              printstate,&
                              gettfinish,&
                              getnpics,&
                              getresolution,&
                              getsilentmode,&
                              getPartNumber,&
                              getHfac,&
                              getGamma,&
                              getresultfile,&
                              getUseDumps,&
                              set_equations,&
                              getState,&
                              settfinish,&
                              getdim,&
                              getLastPrint,&
                              setdtprint, getdtprint
  use iterator,         only: iterate
  use printer,          only: Output, AppendLine
  use errcalc,          only: err_diff_laplace, &
                              err_diff_graddiv, &
                              err_sinxet, &
                              err_shockTube => shockTube, &
                              err_soundwave_d => soundwaveperturbation_density, &
                              err_soundwave_v => soundwaveperturbation_velocity, &
                              err_diff_artvisc => diff_artvisc, &
                              err_alfvenwave => alfvenwave, &
                              err_hcpulse => hcpulse, &
                              err_hcring => hcring, &
                              err_shock12 => diffshock12
  use args,             only: fillargs
  use errteylor,        only: etlaplace => laplace,&
                              etgraddiv => graddiv
  use circuit1,         only: c1destroy => destroy
  use timing,           only: printTimes,&
                              tinit => init, &
                              timedestroy => destroy
  use neighboursearch,  only: getNeibNumbers,&
                              neibDestroy => destroy,&
                              getNeibListL1
  use arrayresize,      only: resize
  use timing,           only: addTime
  use dumper,           only: dump, restore,&
                              dumpclean => clean
  use kernel,           only: initkernel


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
    dt, t, dtout, ltout, tfinish,&
    gamma, hfac, cv = 1., &
    mhdmuzero, difcond, dedt, sumdedt, dedtprev, chi(81)
  character (len=100) :: &
    errfname
  integer :: &
    i, n, iter, s_tt, nusedl1, nusedl2, printlen, silent,&
    ivt, stopiter, resol, difiso, realpartnumb, npics, usedumps,&
    dim, lastnpic

  integer(8) :: &
    tprint

  logical :: thereisdump, needtoswitch

  ! call runcltest()
  ! read*

  print*, "##############################################"
  print*, "#####"
  print*, "#   #   .HELLO WORLD. Setup in progress."
  print*, "#####"
  print*, "##############################################"

  inquire(file='output/fulldump', exist=thereisdump)
  call fillargs()

  t = 0.
  if (thereisdump) then
    call restore('output/fulldump', store, t)
  else
    inquire(file='output/finaldump', exist=thereisdump)
    if (thereisdump) then
      call restore('output/finaldump', store, t)
    else
      call getdim(dim)
      call initkernel(dim)

      call getLastPrint(lastnpic)
      call getnpics(npics)
      call gettfinish(tfinish)
      dtout = (tfinish-t)/(npics - lastnpic)
      call setdtprint(dtout)

      call setupV2(n, cv, store)
    end if
  end if

  call get_equations(s_tt)
  call ginitvar(ivt)
  call getdiffisotropic(difiso)
  call getdiffconductivity(difcond)
  call getmhdmagneticpressure(mhdmuzero)
  call gettfinish(tfinish)
  call getnpics(npics)
  call getsilentmode(silent)
  call getresolution(resol)
  call getGamma(gamma)
  call getHfac(hfac)
  call getresultfile(errfname)
  call getUseDumps(usedumps)
  call getdtprint(dtout)
  call getLastPrint(lastnpic)

  call getdim(dim)
  call initkernel(dim)
  call initBorders()
  call tinit()

  ! select case(s_tt)
  ! case (eeq_hydro, eeq_diffusion, eeq_magnetohydro, eeq_magnetohydrodiffusion, eeq_hydrodiffusion)
  ! ! case (5, 6, 7, 8, 10)
  !   ! 'diff-laplass' ! 'diff-graddiv' ! diff-artvisc
  !   ! call set_stepping(10**dim)
  ! case default
  !   print *, 'Particle turn-off was not set ./src/main.f90:119'
  !   stop
  ! end select

  call printstate()

  allocate(result(100))
  result(:) = 0.
  result(1) = resol
  call getPartNumber(r=realpartnumb)
  result(2) = realpartnumb
  n = realpartnumb
  dt = 0.
  ltout = 0.
  iter = 0.
  allocate(tdu(n))
  allocate(tdh(n))
  allocate(tr(3,n))
  allocate(tv(3,n))
  allocate(ta(3,n))
  allocate(tddt(n))
  allocate(err(1:n))
  allocate(tdb(3,n))
  allocate(sqerr(1:n))

  err(:) = 0.
  stopiter = 0
  sumdedt = 0.

  if (silent == 0) then
    call Output(t, store, sqerr)
    print *, '#  #-------------------------------------------------'
    print *, '#  # Initial setup printed'
    write(*, fmt="(A,I7, A,F12.7, A,F12.7, A,F10.7,A,F10.7,A, A,F12.7, A,F12.7)") &
      " #  #", iter, &
      " | t=", t, &
      " | dt=", dt, &
      " | h=[", minval(store(es_h,1:n)), ":", maxval(store(es_h,1:n)), "]",&
      " | dedt=", dedt*dt,&
      " | S(dedt)=", sumdedt
    print *, '#  #-------------------------------------------------'
  end if

  call checkVarsReady(s_tt, store)

  if (s_tt /= eeq_kd2) then
    call iterate(n, gamma, store, dedt)
  else
    stopiter = 1
  end if
  sumdedt = sumdedt + dedt

  do while(ltout <= t)
    ltout = ltout + dtout
  end do
  do while ((t < tfinish + eps0).and.(stopiter==0))
    ! if (needtoswitch.and.(t > (tfinish/2.))) then
    !   call set_equations('mhd')
    !   s_tt = eeq_magnetohydrodiffusion
    !   needtoswitch = .false.
    !   print*, maxval(store(es_bx:es_bz,1:n)), minval(store(es_bx,1:n))
    ! end if

    select case(s_tt)
    case(eeq_kd2)
      dt = 0.
    case(eeq_hydro)
      dt = .01 * minval(store(es_h,1:n)) / maxval(store(es_c,1:n))
    case(eeq_magnetohydro)
      dt = .1*1e-5*minval(store(es_h,1:n))&
            /maxval(store(es_c,1:n))&
            /maxval(store(es_bx:es_bz,1:n))
      ! dt = .1*minval(store(es_h,1:n))&
      !       /maxval(store(es_c,1:n))
    case(eeq_diffusion, eeq_hydrodiffusion)
      dt = .01 * minval(store(es_den,1:n)) &
              * minval(store(es_c,1:n)) &
              * minval(store(es_h,1:n)) ** 2 &
              / merge(difcond, maxval(store(es_bx:es_bz,1:n)), difiso == 1)
    case (eeq_magnetohydrodiffusion)
      ! dt = .1 * mhd_magneticconstant * minval(den(1:realpartnumb)) * minval(c(1:realpartnumb)) * &
      !           minval(h(1:realpartnumb)) ** 2 / &
      !           merge(maxval(kcf(:,1,1:realpartnumb)), 1., maxval(kcf(:,1,1:realpartnumb))>0)
      ! dt = .3 * minval(store(es_h,1:n)) / maxval(store(es_c,1:n))
      dt = 1.*mhdmuzero&
              * minval(store(es_den,1:n)) &
              * minval(store(es_c,1:n)) &
              * minval(store(es_h,1:n))**2 &
              / difcond
      ! print*, dt
      ! print*, minval(store(es_den,1:n)), minval(store(es_c,1:n)), minval(store(es_p,1:n)), minval(store(es_u,1:n))
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
      if (usedumps == 1) call dump(store, t)
      if (silent == 0) call Output(t, store, err)
      write(*, fmt="(A,I7, A,F12.7, A,F12.7, A,F10.7,A,F10.7,A, A,F12.7, A,F12.7)") &
        " #", iter, &
        " | t=", t, &
        " | dt=", dt, &
        " | h=[", minval(store(es_h,1:n)), ":", maxval(store(es_h,1:n)), "]",&
        " | dedt=", dedt*dt,&
        " | S(dedt)=", sumdedt

      do while(ltout <= t)
        ltout = ltout + dtout
      end do
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
    ! store(es_t,1:n) = store(es_t,1:n) + dt * tddt(:)
    store(es_t,1:n) = store(es_u,1:n)
    store(es_bx:es_bz,1:n) = store(es_bx:es_bz,1:n) + dt * tdb(:,:)
    dedtprev = dedt
    sumdedt = sumdedt + dedt*dt

    call iterate(n, gamma, store, dedt)

    store(es_vx:es_vz,1:n) = &
      store(es_vx:es_vz,1:n) + 0.5*dt*(store(es_ax:es_az,1:n) - ta(:,:))
    store(es_u,1:n) = &
      store(es_u,1:n) + 0.5*dt*(store(es_du,1:n) - tdu(:))
    store(es_h,1:n) = &
      store(es_h,1:n) + 0.5*dt*(store(es_dh,1:n) - tdh(:))
    ! store(es_t,1:n) = &
      ! store(es_t,1:n) + 0.5*dt*(store(es_ddt,1:n) - tddt(:))
    store(es_t,1:n) = store(es_u,1:n)
    store(es_bx:es_bz,1:n) = &
      store(es_bx:es_bz,1:n) + 0.5*dt*(store(es_dbx:es_dbz,1:n) - tdb(:,:))
    sumdedt = sumdedt + 0.5*dt*(dedt - dedtprev)

    t = t + dt
    iter = iter + 1
  end do

  if (s_tt /= eeq_kd2) then
    select case(ivt)
    case (ett_OTvortex, ett_mti, ett_boilingtank, ett_mtilowres)
    case (ett_ring)
      call err_hcring(store, t, err)
    case (ett_pulse)
      call err_hcpulse(store, t, err)
    case (ett_sin3)
      call err_sinxet(store, t, err)
    case (ett_soundwave)
      ! call err_soundwave_v(pos, vel, t, err)
      call err_soundwave_d(store, t, err)
    case (ett_hydroshock)
      call err_shockTube(store, t, err)
    case (ett_alfvenwave)
      call err_alfvenwave(store, t, err)
    case (ett_shock12)
      call err_shock12(store, t, err)
    case default
      print *, ivt, 'Task type was not sen in l2 error evaluation ./src/main.f90:270'
      stop
    end select
    printlen = 5

    call getNeibNumbers(nusedl1, nusedl2)
    result(3) = merge(sqrt(sum(err)/nusedl1), 0., nusedl1 > 0)
  else
    !----------------------------------------!
    !          teylor error evaluation       !
    !----------------------------------------!
    call etlaplace(store, chi)
    result(4) = sum(chi(1:9))/dim
    result(6:14) = chi(1:9)
    printlen = 14
  end if

  if (nusedl1 /= 0) then
    sqerr(:) = sqrt(err(:))
    write(*, fmt="(A,I7, A,F12.7)") &
      " #", iter, &
      " | t=", t
    if (silent == 0) call Output(t, store, sqerr)
    if (usedumps == 1) call dump(store, t)
  end if
  result(5) = hfac
  ! print*, printlen
  call resize(result, printlen, printlen)
  call AppendLine(result, errfname, tprint)

  call printTimes()
  print *, '#####  Results:'
  write(*, "(A, F20.10)") " # #   l2-error: ", result(3)
  write(*, "(A, F20.10)") " # #  chi-error: ", result(4)
  print *, '##############################################'

  call neibDestroy()
  if (s_tt /= eeq_kd2) call c1destroy()
  call timedestroy()
  call dumpclean()
end program

subroutine set_stepping(i)
  use neighboursearch, only: stnb  => setStepsize

  integer, intent(in) :: i

  call stnb(i)
  print *, '# #      step.size:', i
end subroutine
