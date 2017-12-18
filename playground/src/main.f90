program main
  use BC,               only: bcdestroy => destroy
  ! ,&
                              ! ReflectBorders

  use IC,               only: setupIC
  use state,            only: get_tasktype,&
                              ginitvar,&
                              setpartnum,&
                              mhd_magneticconstant
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

  real, allocatable, dimension(:,:,:) :: kcf
  real, allocatable, dimension(:,:)   :: p, v, a, pos, vel, acc, cf, dcf, tcf
  real, allocatable, dimension(:)     :: den, prs, mas, iu, du, om, c, h, dh, &
                                         tdh, err, sqerr, tdu,&
                                         result
  integer, allocatable, dimension(:)  :: ptype, nlista

  real                                :: dt, t, dtout, ltout, tfinish, npic,&
                                         pspc1, gamma,&
                                         sk, chi(81), cv = 1.

  character (len=100) :: itype, errfname, ktype, dtype
  integer             :: n, dim, iter, s_tt, nusedl1, nusedl2, printlen, silent,&
                        ivt, stopiter, resol, i, la

  integer(8)          :: tprint

  ! call runcltest()
  print *, '##############################################'
  print *, '#####'

  call fillargs(dim, pspc1, resol,&
                itype, ktype, dtype, errfname, dtout, npic, tfinish, sk, silent)
  call setupIC(n, sk, gamma, pspc1, resol, pos, vel, acc, &
                mas, den, h, prs, iu, du, cf, kcf, dcf, ptype)

  call setpartnum(n)
  call get_tasktype(s_tt)
  call ginitvar(ivt)

  select case(s_tt)
  case (eeq_hydro, eeq_diffusion, eeq_magnetohydro, eeq_magnetohydrodiffusion)
  case (5, 6, 7, 8, 10)
    ! 'diff-laplass' ! 'diff-graddiv' ! diff-artvisc
    call set_stepping(10**dim)
  case default
    print *, 'Particle turn-off was not set main.f90: line 75.'
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

  call tinit()
  call c1_init(n)
  stopiter = 0
  print *, "Finish time = ", tfinish

  if (silent == 0) then
    print *, 'Initial setup printed'
    call Output(t, ptype, pos, vel, acc, mas, den, h, prs, iu, cf, kcf, sqerr)
  end if

  call checkVarsReady(s_tt, mhd_magneticconstant, den, c, h, prs, iu)

  call iterate(n, sk, gamma, &
              ptype, pos, vel, acc, &
              mas, den, h, dh, om, prs, c, iu, du, &
              cf, dcf, kcf)

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
      dt = .3 * minval(h) / maxval(c)
    case(eeq_magnetohydro)
      dt = .1 * minval(h) / maxval(c)
    case(eeq_diffusion)
      dt = .1 * minval(den) * minval(c) * minval(h) ** 2 / maxval(kcf(:,1,:))
    case (eeq_magnetohydrodiffusion)
      dt = .1 * mhd_magneticconstant * minval(den) * minval(c) * &
                minval(h) ** 2 / merge(maxval(kcf(:,1,:)), 1., maxval(kcf(:,1,:))>0)
      ! dt = .1 * minval(h) / maxval(c)
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
      print *, "#", iter, "t=", t, "dt=", dt, 'min h=', minval(h), 'max h=', maxval(h)
      if ( silent == 0) then
        call Output(t, ptype, pos, vel, acc, mas, den, h, prs, iu, cf, kcf, err)
      end if
      ltout = ltout + dtout
    end if

    if (any(isnan(acc))) then
      print*, 'Acceleration is NAN'
      call Output(t, ptype, pos, vel, acc, mas, den, h, prs, iu, cf, kcf, err)
      exit
    end if

    p(:,:) = pos(:,:)
    v(:,:) = vel(:,:)
    a(:,:) = acc(:,:)
    tdu(:) = du(:)
    tdh(:) = dh(:)
    tcf(:,:) = dcf(:,:)
    kcf(:,3,:) = kcf(:,2,:)
    pos(:,:) = p(:,:)  + dt * v(:,:) + 0.5 * dt * dt * a(:,:)
    ! call ReflectBorders()
    vel(:,:) = v(:,:)  + dt * a(:,:)
    iu(:)    = iu(:)   + dt * du(:)
    h(:)     = h(:)    + dt * dh(:)
    cf(:,:)  = cf(:,:) + dt * dcf(:,:)
    kcf(:,1,:) = kcf(:,1,:) + dt * kcf(:,2,:)
    call iterate(n, sk, gamma, &
                ptype, pos, vel, acc, &
                mas, den, h, dh, om, prs, c, iu, du, &
                cf, dcf, kcf)
    vel(:,:) = vel(:,:) + 0.5 * dt * (acc(:,:)  - a(:,:))
    iu(:)    = iu(:)    + 0.5 * dt * (du(:)     - tdu(:))
    h(:)     = h(:)     + 0.5 * dt * (dh(:)     - tdh(:))
    cf(:,:)  = cf(:,:)  + 0.5 * dt * (dcf(:,:)  - tcf(:,:))
    kcf(:,1,:) = kcf(:,1,:) + 0.5 * dt * (kcf(:,2,:)  - kcf(:,3,:))
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
  case (-1, 2, 3, 4, 5, ett_OTvortex)
    ! ! diff-laplace ! chi-laplace
    ! call etlaplace(pos, mas, den, h, chi)
    ! result(4) = sum(chi(1:9))/dim
    ! result(6:14) = chi(1:9)
    ! printlen = 14

    ! ! diff-graddiv ! chi-graddiv
    ! call etgraddiv(pos, mas, den, h, chi)
    ! result(4) = sum(chi)/dim/(2*dim - 1)
    ! result(6:86) = chi(1:81)
    ! printlen = 86
    printlen = 5
  case (ett_sin3)
    call err_sinxet(pos, cf, t, err)
    printlen = 5
  case (ett_soundwave)
    ! call err_soundwave_v(pos, vel, t, err)
    call err_soundwave_d(pos, den, t, err)
    printlen = 5
  case (ett_hydroshock)
    call err_shockTube(ptype, pos, den, t, err)
    printlen = 5
  case (ett_alfvenwave)
    call err_alfvenwave(pos, cf, t, err)
    printlen = 5
  case default
    print *, ivt, 'Task type was not sen in l2 error evaluation main.f90: line 229.'
    stop
  end select

  call getNeibNumbers(nusedl1, nusedl2)
  result(3) = merge(sqrt(sum(err)/nusedl1), 0., nusedl1 > 0)
  !----------------------------------------!
  !          teylor error evaluation       !
  !----------------------------------------!
  if (nusedl1 /= 0) then
    sqerr(:) = sqrt(err(:))
    print *, iter, t, 0.
    if (silent == 0) then
      call Output(t, ptype, pos, vel, acc, mas, den, h, prs, iu, cf, kcf, sqerr)
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
