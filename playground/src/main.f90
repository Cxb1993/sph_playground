program main
  ! use cudafor
  use BC,               only: initBorders
  use IC,               only: setupV2
  use placer,           only: place
  use state,            only: clearState,&
                              get_equations,&
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
                              err_alfvenwave => alfvenwave,&
                              err_hcpulse => hcpulse,&
                              err_hcring => hcring,&
                              err_hcshock12 => hcshock12
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
  use arrayresize,      only: resize!, realocnum
  use timing,           only: addTime
  use dumper,           only: dump, restore,&
                              dumpclean => clean
  use kernel,           only: initkernel,&
                              getcndiff,&
                              getcnhydro
  use errprinter,       only: warning, error

  use timestep,         only: calctimestep => calc

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
    dedt, sumdedt, dedtprev, chi(81),&
    avdt
  character (len=100) :: &
    errfname
  integer :: &
    n, iter, s_tt, nusedl1, nusedl2, printlen, silent,&
    ivt, stopiter, resol, realpartnumb, npics, usedumps,&
    dim, lastnpic

  integer(8) :: &
    tprint

  logical :: thereisdump

  ! call runcltest()
  ! read*

  call clearState()

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
      if (dtout <= 0.) then
        dtout = 1.
        call warning("Dirty hack was used here", "dtout", __FILE__, __LINE__)
      end if
      call setdtprint(dtout)

      ! hopefully that is just for some time to supprt depricated ivt methods
      call ginitvar(ivt)
      if (ivt == e_none) then
        call place(store)
      else
        call setupV2(n, cv, store)
      end if
    end if
  end if

  call get_equations(s_tt)
  call ginitvar(ivt)
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
  ! call diffinit()

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
  dedt = 0.

  if (usedumps == 1) call dump(store, t)
  if (silent == 0) call Output(t, store, sqerr)
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

  if (tfinish > 0.) then
    call checkVarsReady(s_tt, store)
    if (s_tt /= eeq_kd2) then
      call iterate(n, gamma, store, dedt)
    else
      stopiter = 1
    end if
  else
    stopiter = 1
  end if

  sumdedt = sumdedt + dedt

  do while(ltout <= t)
    ltout = ltout + dtout
  end do

  do while ((t < tfinish + eps0).and.(stopiter==0))
    call calctimestep(store, dt)

    if (iter == 0) then
      avdt = dt
    else
      avdt = avdt + dt
    end if
    if (t + dt > tfinish) then
      dt = tfinish - t
      stopiter = 1
    end if

    if (t >= ltout) then
      if (usedumps == 1) call dump(store, t)
      if (silent == 0) call Output(t, store, err)
      write(*, fmt="(A,I7, A,F12.7, A,F12.7, A,F10.7,A,F10.7,A, A,F12.7, A,F12.7)") &
        "#", iter, &
        " | t=", t, &
        " | dt=", merge(avdt/iter, avdt, iter>0), &
        " | h=[", minval(store(es_h,1:n)), ":", maxval(store(es_h,1:n)), "]",&
        " | dedt=", dedt*dt,&
        " | S(dedt)=", sumdedt

      do while(ltout <= t)
        ltout = ltout + dtout
      end do
    end if
    if (t + dt > ltout) then
      dt = ltout - t
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

    store(es_vx:es_vz,1:n) = (tv(:,:)  + dt * ta(:,:))/(1.03)
    store(es_u,1:n) = store(es_u,1:n) + dt * tdu(:)
    store(es_h,1:n) = store(es_h,1:n) + dt * tdh(:)
    store(es_t,1:n) = store(es_u,1:n)
    store(es_bx:es_bz,1:n) = store(es_bx:es_bz,1:n) + dt * tdb(:,:)
    dedtprev = dedt
    sumdedt = sumdedt + dedt*dt

    call iterate(n, gamma, store, dedt)

    store(es_vx:es_vz,1:n) = &
      (store(es_vx:es_vz,1:n) + 0.5*dt*(store(es_ax:es_az,1:n) - ta(:,:)))/(1.03)
    store(es_u,1:n) = &
      store(es_u,1:n) + 0.5*dt*(store(es_du,1:n) - tdu(:))
    store(es_h,1:n) = &
      store(es_h,1:n) + 0.5*dt*(store(es_dh,1:n) - tdh(:))
    store(es_bx:es_bz,1:n) = &
      store(es_bx:es_bz,1:n) + 0.5*dt*(store(es_dbx:es_dbz,1:n) - tdb(:,:))
    sumdedt = sumdedt + 0.5*dt*(dedt - dedtprev)
    store(es_t,1:n) = store(es_u,1:n)

    t = t + dt
    iter = iter + 1
  end do

  if (tfinish > 0.) then
    if (s_tt /= eeq_kd2) then
      select case(ivt)
      case (ett_OTvortex, ett_mti, ett_boilingtank, ett_mtilowres)
      case (ett_shock12)
        call err_hcshock12(store, t, err)
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
      case (e_none)
        call warning('L2 error evaluation: none was used as a initial problem type', '', __FILE__, __LINE__)
      case default
        call error('Task type was not sen in L2 error evaluation', '', __FILE__, __LINE__)
      end select
      printlen = 5

      call getNeibNumbers(nusedl1, nusedl2)
      result(3) = merge(sum(err)/nusedl1, 0., nusedl1 > 0)
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

  end if

  call printTimes()
  if (tfinish > 0) then
    print *, '#####  Results:'
    write(*, "(A, F20.10)") " # #   L1-error: ", result(3)
    write(*, "(A, F20.10)") " # #  chi-error: ", result(4)
    print *, '##############################################'
  end if
  ! call neibDestroy()
  ! if (s_tt /= eeq_kd2) call c1destroy()
  ! call timedestroy()
  call dumpclean()
end program

subroutine set_stepping(i)
  use neighboursearch, only: stnb  => setStepsize

  integer, intent(in) :: i

  call stnb(i)
  print *, '# #      step.size:', i
end subroutine
