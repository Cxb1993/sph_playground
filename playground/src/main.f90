program main
  ! use cudafor
  use omp_lib
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
                              getUseDumps,&
                              set_equations,&
                              getState,&
                              settfinish,&
                              getdim,&
                              getLastPrint,&
                              setdtprint, getdtprint,&
                              getStateVal,&
                              setStateVal,&
                              getEqComponent
  use iterator,         only: iterate
  use printer,          only: AppendLine, handleOutput
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

  use rad_exchange,     only:update_radenergy

  use errprinter,       only: warning, error

  use timestep,         only: calctimestep => calc

  use sts_integrator,   only: sts_test => test,&
                              sts_integrate,&
                              sts_init => init

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
    dt, t, dtout, dtoutaim,ltout, tfinish,&
    gamma, hfac, cv = 1., &
    dedt, sumdedt, dedtprev, chi(81),&
    sumdt, phystime, phystimeprev, physdt,&
    dthydro, dtdiff, dtfld,&
    sts_s_sum,sts_s_iter,&
    sts_dt_min,sts_classic,sts_classic_sum,sts_classic_iter
  integer :: &
    i, n, iter, s_tt, nusedl1, nusedl2, printlen, silent,&
    ivt, stopiter, resol, realpartnumb, fixedpartnumb, npics, usedumps,&
    dim,lastnpic,eqSet(eqs_total),sts_s,sts_type,sts_fixeds
  integer(8) :: &
    tprint
  logical :: thereisdump, do_sts=.false.
  ! call runcltest()
  ! read*

  call clearState()

  print*, "##############################################"
  print*, "#####"
  print*, "#   #   .HELLO WORLD. Setup is in progress."
  print*, "#####"
  print*, "##############################################"

  t = 0.
  inquire(file='output/dump_full.h5', exist=thereisdump)
  if (thereisdump) then
    call restore(store)
    call getStateVal(ec_time, t)
  else
    call fillargs()

    call getStateVal(ec_dim, dim)
    call initkernel(dim)
    call getStateVal(ec_lastprint, lastnpic)
    call getStateVal(ec_npics, npics)
    call getStateVal(ec_tfinish, tfinish)

    dtout = (tfinish-t)/(npics - lastnpic)
    if (dtout <= 0.) then
      dtout = 1.
      call warning("Dirty hack was used here", "dtout", __FILE__, __LINE__)
      call warning("And I even don't remember why enymore", "dtout", __FILE__, __LINE__)
    end if
    call setStateVal(ec_dtprint, dtout)

    ! hopefully that is just for some time to supprt depricated ivt methods
    call getStateVal(ec_ics, ivt)
    if (ivt == e_none) then
      call place(store)
    else
      call setupV2(n, cv, store)
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
  call getStateVal(ec_fixedpn, fixedpartnumb)
  call getStateVal(ec_realpn, realpartnumb)
  result(2) = realpartnumb

  call getEqComponent(eqSet)

  n = realpartnumb
  do i = realpartnumb,size(store,dim=2)
    if (store(es_type,i) == ept_fixedreal) n = n + 1
  end do
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

  if (eqSet(eqs_sts) == 1) then
    call sts_init(n)
  endif

  err(:) = 0.
  stopiter = 0
  sumdedt = 0.
  dedt = 0.
  sumdt = 0.
  physdt = 0.
  print *, '#  #-------------------------------------------------'
  print *, '#  # Initial Setup'
  call handleOutput(iter, n, t, sumdt, dedt, dt, sumdedt, physdt, store, err)
  print *, '#  #-------------------------------------------------'
#ifdef _OPENMP
  phystimeprev = omp_get_wtime()
#endif

  if (tfinish > 0.) then
    call checkVarsReady(s_tt, store)
    if (s_tt /= eeq_kd2) then
      ! print*, 1
      call iterate(n, store, dedt)
    else
      stopiter = 1
    end if
  else
    stopiter = 1
  end if

  sumdedt = sumdedt + dedt

  do while(ltout < t)
    ltout = ltout + dtout
  end do
  dtoutaim = dtout

  sts_s_sum = 0.
  sts_s_iter = 0.
  sts_classic_sum = 0.
  sts_classic_iter = 0.
  call getStateVal(ec_ststype, sts_type)
  call getStateVal(ec_stsfixeds, sts_fixeds)

  do while ((t < tfinish + eps0).and.(stopiter==0))
    call calctimestep(store, dthydro, dtdiff, dtfld, dt)
    if (eqSet(eqs_sts) == 1) then
      do_sts = .false.
      ! >>>>> STS force hydro (check unstable)
      ! do_sts = .false.
      ! dt = dthydro
      sts_dt_min = min(dtdiff,dtfld)
      if (sts_type == ests_fixeds) then
        ! >>>>> STS fixed stages, varied timestep
        sts_s = sts_fixeds
        dthydro = (sts_s*sts_s + sts_s - 2.)*0.25*sts_dt_min
        sts_classic = dthydro/sts_dt_min
        dt = dthydro
        do_sts = .true.
      else if (sts_type == ests_auto) then
        ! >>>>> STS auto stages choice, hydro timestep
        if (dthydro > 0.) then
          sts_classic = dthydro/sts_dt_min
          sts_s = int((-1.+sqrt(1.-4.*1.*(-(sts_classic*4.-2.))))*0.5)+1
          if ((sts_classic > 1.).and.(2*sts_s <= sts_classic)) then
            dt = dthydro
            do_sts = .true.
          endif
        endif
      endif
      ! >>>>> STS status bar
      sts_classic_sum  = sts_classic_sum  + max(1.,sts_classic)
      sts_classic_iter = sts_classic_iter + 1
      sts_s_sum  = sts_s_sum  + sts_s
      sts_s_iter = sts_s_iter + 1
    endif

    if (t + dt > tfinish) then
      dt = tfinish - t
      stopiter = 1
    end if
    sumdt = sumdt + dt

    if (t >= ltout) then
#ifdef _OPENMP
      phystime = omp_get_wtime()
#endif
      physdt = phystime - phystimeprev
      call handleOutput(iter, n, t, sumdt, dedt, dt, sumdedt, physdt, store, err)
      if (eqSet(eqs_sts) == 1) then
        if (do_sts) then
          print*, "#        | V | STS-ON:  sts steps = ",int(sts_s_sum/sts_s_iter),&
            dthydro, &
            ": classic steps = ",&
            int(sts_classic_sum/sts_classic_iter),sts_dt_min
        else
          print*, "#        | X | STS-OFF: sts steps = ",int(sts_s_sum/sts_s_iter),&
          dthydro, &
          ": classic steps = ",&
          int(sts_classic_sum/sts_classic_iter),sts_dt_min
        endif
        sts_s_sum = 0.
        sts_s_iter = 0.
        sts_classic_sum = 0.
        sts_classic_iter = 0.
      endif
      phystimeprev = phystime
      do while(ltout <= t)
        ltout = ltout + dtout
      end do
    end if
    if (t + dt > ltout) then
      if (dtout < dt) then
        ltout = ltout - dtout + dt
      else
        dt = ltout - t
      endif
    end if

    ! if (do_sts) then
    !   call sts_integrate(n,sts_s,store,dt/2.)
    ! endif
    ! if (eqSet(eqs_radexch)==1) then
    !   call update_radenergy(n,store,dt/2.)
    ! endif

    tr(:,:)  = store(es_rx:es_rz,1:n)
    tv(:,:)  = store(es_vx:es_vz,1:n)
    ta(:,:)  = store(es_ax:es_az,1:n)
    tdu(:)   = store(es_du,1:n)
    tdh(:)   = store(es_dh,1:n)
    if (.not.do_sts) then
      tddt(:)  = store(es_ddt,1:n)
    endif
    tdb(:,:) = store(es_dbx:es_dbz,1:n)

    store(es_rx:es_rz,1:n) = store(es_rx:es_rz,1:n)  + &
        dt * tv(:,:) + 0.5 * dt * dt * ta(:,:)

    store(es_vx:es_vz,1:n) = tv(:,:)  + dt * ta(:,:)
    store(es_u,1:n) = store(es_u,1:n) + dt * tdu(:)
    store(es_h,1:n) = store(es_h,1:n) + dt * tdh(:)
    if (.not.do_sts) then
      store(es_t,1:n) = store(es_t,1:n) + dt * tddt(:)
    endif
    store(es_bx:es_bz,1:n) = store(es_bx:es_bz,1:n) + dt * tdb(:,:)
    dedtprev = dedt
    sumdedt = sumdedt + dedt*dt

    call iterate(n, store, dedt)

    store(es_vx:es_vz,1:n) = &
      store(es_vx:es_vz,1:n) + 0.5*dt*(store(es_ax:es_az,1:n) - ta(:,:))
    ! store(es_vx:es_vz,1:n) = &
    !   (store(es_vx:es_vz,1:n) + 0.5*dt*(store(es_ax:es_az,1:n) - ta(:,:)))/(1.03)
    store(es_u,1:n) = &
      store(es_u,1:n) + 0.5*dt*(store(es_du,1:n) - tdu(:))
    store(es_h,1:n) = &
      store(es_h,1:n) + 0.5*dt*(store(es_dh,1:n) - tdh(:))
    store(es_bx:es_bz,1:n) = &
      store(es_bx:es_bz,1:n) + 0.5*dt*(store(es_dbx:es_dbz,1:n) - tdb(:,:))
    sumdedt = sumdedt + 0.5*dt*(dedt - dedtprev)
    if (.not.do_sts) then
      store(es_t,1:n) = &
        store(es_t,1:n) + 0.5*dt*(store(es_ddt,1:n) - tddt(:))
    endif

    if (do_sts) then
      call sts_integrate(n,sts_s,store,dt)
    endif
    if (eqSet(eqs_radexch)==1) then
      call update_radenergy(n,store,dt)
    endif

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
      write(*, fmt="(A,I7, A,ES10.4)") &
        " #", iter, &
        " | t=", t
#ifdef _OPENMP
        phystime = omp_get_wtime()
#endif
        physdt = phystime - phystimeprev
        call handleOutput(iter, n, t, sumdt, dedt, dt, sumdedt, physdt, store, err)
        phystimeprev = phystime
    end if
    result(5) = hfac
    ! print*, printlen
    call resize(result, printlen, printlen)
    call AppendLine(result, 'result.info', tprint)
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
