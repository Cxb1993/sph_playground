module circuit1
  use const
  use omp_lib
  use timing,           only: addTime
  use state,            only: getdim, &
                              getddwtype, &
                              gcoordsys, &
                              ginitvar, &
                              getAdvancedDensity,&
                              getPartNumber,&
                              gethfac,&
                              getEqComponent
  use kernel,           only: get_krad, &
                              get_dw_dh, &
                              nw, &
                              get_w
  use neighboursearch,  only: getneighbours, &
                              getNeibListL1, &
                              getNeibListL2, &
                              xgetneighnumber, &
                              xgetneighindex
  use BC,               only: getCrossRef

  implicit none

  public :: c1, destroy

  private
  save
    real, allocatable :: &
      hk(:), hkp1(:), eps(:), dk(:), dkp1(:), ok(:), okp1(:)
    integer(8) :: start=0, finish=0
    integer :: &
      initialised=0,&
      realpartnumb, fixedpartnumb, totreal, eqSet(eqs_total), ad

contains
  subroutine init()
    call getPartNumber(r=realpartnumb, f=fixedpartnumb)
    totreal = realpartnumb + fixedpartnumb
    allocate(hk(totreal))
    allocate(hkp1(totreal))
    allocate(dk(totreal))
    allocate(dkp1(totreal))
    allocate(ok(totreal))
    allocate(okp1(totreal))
    allocate(eps(totreal))
    call getEqComponent(eqSet)
    call getAdvancedDensity(ad)
    initialised = 1
  end subroutine init

  subroutine c1(store)
    real, allocatable, intent(inout) :: store(:,:)

    if (initialised == 0) call init()

    if (ad == 1) then
      call c1advanced(store)
    else
      call c1simple(store)
      store(es_om,:) = 1.
    end if
  end subroutine

  subroutine destroy()
    deallocate(hk)
    deallocate(hkp1)
    deallocate(dk)
    deallocate(dkp1)
    deallocate(ok)
    deallocate(okp1)
    deallocate(eps)
  end subroutine

  subroutine c1advanced(store)
    real, allocatable, intent(inout) :: store(:,:)

    real :: &
      w, dwdh, r(3), dr, r2, dfdh, fh, hn, &
      allowerror, maxinterr, currinterr, &
      mb, ma, ha0, ha, hb, da, db, ra(3), rb(3), oma, omb, &
      dtadx(3), hfac, &
      nwa(3), nwb(3), ta, tb
    integer :: &
      i, j, rj, la, lb, dim, iter, ktp
    integer(8) :: &
      t0, tneib
    integer, allocatable :: &
      nlista(:), nlistb(:)

    call system_clock(start)

    call getdim(dim)
    call getddwtype(ktp)
    call gethfac(hfac)

    call getNeibListL1(nlista)
    ! print*, size(nlista)
    ! call getNeibListL2(nlista)
    ! print*, size(nlista)
    ! read*

    allowerror = 1e-8
    hkp1(:) = store(es_h,  1:totreal)
    dkp1(:) = store(es_den,1:totreal)
    okp1(:) = store(es_om, 1:totreal)
    eps(:)  = 0.
    eps(nlista) = 1.
    iter = 0
    tneib = 0.
    maxinterr  = 0.
    currinterr = 0.
    do while ((maxval(eps(:), mask=(eps>0)) > allowerror) .and. (iter < 100))
      maxinterr  = 0.
      iter = iter + 1
      hk(:) = hkp1(:)
      dk(:) = dkp1(:)
      ok(:) = okp1(:)
      !$omp parallel do default(none)&
      !$omp private(r, dr, dwdh, w, dfdh, fh, hn, j, rj, i, la, lb, r2, t0, nlistb)&
      !$omp private(currinterr, mb, ma, da, db, ra, rb, dtadx)&
      !$omp private(nwa, nwb, ta, tb, ha0, ha, hb, oma, omb)&
      !$omp shared(eps, allowerror, dim, hfac, ktp, maxinterr)&
      !$omp shared(store, nlista, dk, dkp1, hk, hkp1, ok, okp1)&
      !$omp shared(nw, eqSet)&
      !$omp reduction(+:tneib)
      do la = 1, size(nlista)
        i = nlista(la)
        if (eps(i) > allowerror) then
          currinterr = 0.
          ma  = store(es_m,i)
          ha0 = store(es_h,i)
          ra(:) = store(es_rx:es_rz,i)
          ta  = store(es_t,i)
          ha  = hk(i)
          da  = dk(i)
          oma = ok(i)
          dtadx(:) = 0.
          dkp1(i)  = 0.
          hkp1(i)  = 0.
          okp1(i)  = 0.

          tneib = tneib + t0
          ! print*, i, xgetneighnumber(i)
          ! do lb = 1, xgetneighnumber(i)
          !   print*, xgetneighindex(i, lb)
          ! end do
          ! read
          do lb = 1, xgetneighnumber(i)
            j = xgetneighindex(i, lb)
          ! call getneighbours(i, nlistb, t0)
          ! do lb = 1, size(nlistb)
          !   j = nlistb(lb)

            rj = getCrossRef(j)
            rb(:) = store(es_rx:es_rz, j)
            r(:)  = ra(:) - rb(:)
            r2 = dot_product(r(:),r(:))
            dr = sqrt(r2)
            mb = store(es_m,rj)
            tb = store(es_t,rj)
            db = dk(rj)
            hb = hk(rj)
            omb = ok(rj)

            call get_dw_dh(dr, ha, dwdh)
            call get_w(dr, ha, w)
            dkp1(i) = dkp1(i) + mb * w
            okp1(i) = okp1(i) + mb * dwdh

            if (eqSet(eqs_diff) == 1) then
              if (ktp == esd_2nw_ds) then
                call nw(r(:), ra(:), rb(:), dr, ha, nwa)
                dtadx(:) = dtadx(:) + mb/da/oma*(tb - ta)*nwa(:)
              else if (ktp == esd_2nw_sd) then
                ! print*, nlistb
                call nw(r(:), ra(:), rb(:), dr, ha, nwa)
                call nw(r(:), ra(:), rb(:), dr, hb, nwb)
                dtadx(:) = dtadx(:) + mb*(ta/oma/da/da*nwa(:) + tb/omb/db/db*nwb(:))
              end if
            end if
            ! if (i == 1) then
            !   print*, i, j
            !   print*, da, db, dkp1(i)
            !   read*
            ! end if
          end do
          ! ---------------------------------------------------------!
          !      There is no particle itself in neighbour list       !
          ! ---------------------------------------------------------!
          call get_w(0., ha, w)
          call get_dw_dh(0., ha, dwdh)
          dkp1(i) = dkp1(i) + ma*w
          okp1(i) = okp1(i) + ma*dwdh
          ! that is right according to Price2008, not phantom paper
          okp1(i) = 1. - okp1(i)*(-ha/(dim*dkp1(i)))
          dfdh = -dim*dkp1(i)*okp1(i)/ha
          ! print*, 3
          fh   = ma*(hfac/ha)**dim - dkp1(i)
          ! print*, 4
          hn   = ha - fh/dfdh
          if (hn <= 0.) then
            error stop "Negative smoothing length in c1 advanced"
          end if
          eps(i) = abs(hn - ha)/ha0
          hkp1(i) = hn
          store(es_dtdx:es_dtdz,i) = dtadx(:)
        end if
      end do
      !$omp end parallel do
    end do
    ! print*, iter
    store(es_h,  1:realpartnumb) = hkp1(1:realpartnumb)
    store(es_den,1:realpartnumb) = dkp1(1:realpartnumb)
    store(es_om, 1:realpartnumb) = okp1(1:realpartnumb)

    if (iter > 10) then
      print*, "Warn: density NR: solution took ", iter, "iterations, with max norm error", maxval(eps, mask=(eps>0))
    end if

    call system_clock(finish)
    call addTime(' circuit1', finish - start - tneib)
  end subroutine

! Direct density summation
  subroutine c1simple(store)
    real, allocatable, intent(inout) :: store(:,:)
    real :: &
      w, r(3), dr, ra(3), rb(3), ha, hb, ma, mb, da, db, hfac,&
      dtadx(3), nwa(3), nwb(3), ta, tb
    integer :: &
      i, j, rj, la, lb, dim, ktp
    integer, allocatable             :: nlista(:), nlistb(:)
    integer(8)                       :: t0, tneib

    call system_clock(start)
    call getNeibListL1(nlista)
    call getdim(dim)
    call gethfac(hfac)
    call getddwtype(ktp)

    dk(:) = store(es_den,1:totreal)
    hk(:) = store(es_h,1:totreal)

    tneib = 0.
    !$omp parallel do default(none)&
    !$omp private(ra, rb, ha, hb, ma, mb, da, db)&
    !$omp private(r, dr, w, j, rj, i, la, lb, nlistb, t0)&
    !$omp private(dtadx, nwa, nwb, ta, tb)&
    !$omp shared(store, hfac, dim,nlista, dk, hk, nw, eqSet, ktp)&
    !$omp reduction(+:tneib)
    do la = 1, size(nlista)
      i = nlista(la)

      dk(i) = 0.
      hk(i) = 0.
      dtadx(:) = 0.
      ra(:) = store(es_rx:es_rz,i)
      ha = store(es_h,i)
      ma = store(es_m,i)
      da = store(es_den,i)
      ta = store(es_t,i)
      ! call getneighbours(i, nlistb, t0)
      ! tneib = tneib + t0
      ! do lb = 1, size(nlistb)
      !   j = nlistb(lb)
      do lb = 1, xgetneighnumber(i)
        j = xgetneighindex(i, lb)
        rj = getCrossRef(j)
        rb(:) = store(es_rx:es_rz, j)
        hb = store(es_h,i)
        mb = store(es_m, rj)
        db = store(es_den, rj)
        tb = store(es_t,rj)

        r(:) = ra(:)-rb(:)
        dr = sqrt(dot_product(r(:),r(:)))
        call get_w(dr, ha, w)
        dk(i) = dk(i) + mb * w
        if (eqSet(eqs_diff) == 1) then
          if (ktp == esd_2nw_ds) then
            call nw(r(:), ra(:), rb(:), dr, ha, nwa)
            dtadx(:) = dtadx(:) + mb/da*(tb - ta)*nwa(:)
          else if (ktp == esd_2nw_sd) then
            ! print*, nlistb
            call nw(r(:), ra(:), rb(:), dr, ha, nwa)
            call nw(r(:), ra(:), rb(:), dr, hb, nwb)
            dtadx(:) = dtadx(:) + mb*(ta/da/da*nwa(:) + tb/db/db*nwb(:))
          end if
        end if
      end do
      call get_w(0., ha, w)
      dk(i) = dk(i) + ma * w
      hk(i) = hfac * (ma / dk(i))**(1./dim)
      store(es_dtdx:es_dtdz,i) = dtadx(:)
    end do
    !$omp end parallel do
    store(es_den,1:realpartnumb) = dk(1:realpartnumb)
    store(es_h,1:realpartnumb)   = hk(1:realpartnumb)

    call system_clock(finish)
    call addTime(' circuit1', finish - start - tneib)
  end subroutine
end module
