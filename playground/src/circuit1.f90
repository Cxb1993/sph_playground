module circuit1
  use const
  use omp_lib
  use timing,           only: addTime
  use state,            only: getdim, &
                              getddwtype, &
                              gcoordsys, &
                              ginitvar, &
                              getAdvancedDensity,&
                              getPartNumber
  use kernel,           only: get_krad, &
                              get_dw_dh, &
                              nw, &
                              get_w
  use neighboursearch,  only: getneighbours,&
                              getNeibListL1,&
                              getNeibListL2
  use BC,               only: getCrossRef

  implicit none

  public :: c1_init, c1, destroy

  private
  save
    real, allocatable :: slnint(:), resid(:), dennew(:)
    integer(8) :: start=0, finish=0
    integer :: realpartnumb

contains
  subroutine c1_init()
    call getPartNumber(r=realpartnumb)
    allocate(slnint(realpartnumb))
    allocate(resid(realpartnumb))
    allocate(dennew(realpartnumb))
  end subroutine c1_init

  subroutine c1(store, hfac)
    real, allocatable, intent(inout) :: store(:,:)
    real, intent(in) :: hfac

    integer :: ad

    call getAdvancedDensity(ad)

    if (ad == 1) then
      call c1advanced(store, hfac)
    else
      call c1simple(store, hfac)
      store(es_om,:) = 1.
    end if
  end subroutine

  subroutine destroy()
    deallocate(slnint)
    deallocate(resid)
    deallocate(dennew)
  end subroutine

  subroutine c1advanced(store, hfac)
    real, allocatable, intent(inout) :: store(:,:)
    real, intent(in)     :: hfac

    real :: &
      w, dwdh, r(3), dr, r2, dfdh, fh, hn, nwa(3), &
      allowerror, maxinterr, currinterr, &
      mb, ma, ha, da, db, ra(3), rb(3), oma, dtadx(3)
    integer :: &
      i, j, rj, la, lb, dim, iter, ktp
    integer(8)           :: t0, tneib
    integer, allocatable :: nlista(:), nlistb(:)

    call system_clock(start)

    call getdim(dim)
    call getddwtype(ktp)

    call getNeibListL1(nlista)
    allowerror = 1e-8
    slnint(:) = store(es_h,1:realpartnumb)
    dennew(:) = store(es_den,1:realpartnumb)
    resid(:)  = 0.
    resid(nlista) = 1.
    iter = 0
    tneib = 0.
    maxinterr  = 0.
    currinterr = 0.

    do while ((maxval(resid, mask=(resid>0)) > allowerror) .and. (iter < 100))
      maxinterr  = 0.
      iter = iter + 1
      !$omp parallel do default(none)&
      !$omp private(r, dr, dwdh, w, dfdh, fh, hn, j, rj, i, la, lb, r2, t0, nlistb)&
      !$omp private(nwa, currinterr, mb, ma, ha, da, db, ra, rb, oma, dtadx)&
      !$omp shared(resid, allowerror, dim, hfac, ktp, maxinterr)&
      !$omp shared(store, nlista, dennew, slnint)&
      !$omp shared(nw)&
      !$omp reduction(+:tneib)
      do la = 1, size(nlista)
        i = nlista(la)
        if (resid(i) > allowerror) then
          dennew(i) = 0.
          oma = 0.
          currinterr = 0.
          dtadx(:) = 0.

          ma = store(es_m,i)
          ha = slnint(i)
          ra(:) = store(es_rx:es_rz,i)
          da = store(es_den,i)
          call getneighbours(i, nlistb, t0)
          tneib = tneib + t0
          ! print*, nlista
          ! print*, "a=", i, "   neibs=", nlistb, pos(:,i)
          ! read*
          do lb = 1, size(nlistb)
            j = nlistb(lb)
            rj = getCrossRef(j)

            rb(:) = store(es_rx:es_rz, j)
            r(:) = ra(:) - rb(:)
            r2 = dot_product(r(:),r(:))
            ! print*, -3
            dr = sqrt(r2)
            mb = store(es_m, rj)
            db = store(es_den, rj)
            ! print*, -2
            call get_dw_dh(dr, ha, dwdh)
            ! print*, -1
            call get_w(dr, ha, w)
            ! print*,  0
            dennew(i) = dennew(i) + mb * w
            oma = oma + mb * dwdh
            currinterr = currinterr + mb/db * w
            if (ktp == esd_2nw) then
              call nw(r(:), ra(:), rb(:), dr, ha, nwa)
              dtadx(:) = dtadx(:) + mb/db*(store(es_t, rj) - store(es_t, i))*nwa(:)
            end if
          end do
          ! ---------------------------------------------------------!
          !      There is no particle itself in neighbour list       !
          ! ---------------------------------------------------------!
          call get_w(0., ha, w)
          dennew(i) = dennew(i) + ma * w
          currinterr = currinterr + ma/da * w
          if (currinterr > maxinterr) then
            maxinterr = currinterr
          end if

          call get_dw_dh(0., ha, dwdh)
          oma = oma + ma * dwdh
          ! -(**)----------------------------------------------------!
          ! print*,'c1', 4, om(i), mas(i), dwdh
          oma = 1. - oma * (-ha / (dim * dennew(i)))
          ! print*,'c1', 5, i, om(i), slnint(i), dim, den(i)
          dfdh = - dim * dennew(i) * oma / ha
          ! print*,'c1', 7, den(i), om(i), slnint(i)
          fh  = ma * (hfac / ha) ** dim - dennew(i)
          ! print*,'c1', 8, dfdh
          hn = ha - fh / dfdh
          if (hn <= 0.) then
            error stop "Negative smoothing length in c1 advanced"
          end if
          ! print*,'c1', 9
          resid(i) = abs(hn - ha) / store(es_h,i)
          slnint(i) = hn
          store(es_om,i) = oma
          store(es_dtdx:es_dtdz, i) = dtadx(:)
        end if
        ! print *,'c1', dennew(i), slnint(i), om(i)
        ! read*
      end do
      !$omp end parallel do
      store(es_den,1:realpartnumb) = dennew(:)
    end do
    ! print*, '2'
    if (iter > 10) then
      print*, "Warn: density NR: solution took ", iter, "iterations, with max norm error", maxval(resid, mask=(resid>0))
    end if

    ! if (abs(1. - maxinterr) > 0.1) then
    !   print*, "Warn: density NR: kernel integral condition does not fulfilled Int(V_{ab}*W_{ab}) = ", maxinterr
    ! end if
    ! read*
    store(es_h,1:realpartnumb) = slnint(:)
    call system_clock(finish)
    call addTime(' circuit1', finish - start - tneib)
  end subroutine

! Direct density summation
  subroutine c1simple(store, hfac)
    real, allocatable, intent(inout) :: store(:,:)
    real,              intent(in)    :: hfac
    real :: &
      w, r(3), dr, currinterr, ra(3), rb(3), ha, ma, mb, da, db
    integer :: &
      i, j, la, lb, dim
    integer, allocatable             :: nlista(:), nlistb(:)
    integer(8)                       :: t0, tneib

    call system_clock(start)
    call getNeibListL2(nlista)
    call getdim(dim)

    realpartnumb = size(nlista)
    dennew(:) = store(es_den,1:realpartnumb)
    currinterr = 0.

    tneib = 0.
    !$omp parallel do default(none)&
    !$omp private(ra, rb, ha, ma, mb, da, db)&
    !$omp private(r, dr, w, j, i, la, lb, nlistb, t0, currinterr)&
    !$omp shared(store, hfac, dim,nlista, dennew)&
    !$omp reduction(+:tneib)
    do la = 1, size(nlista)
      i = nlista(la)

      dennew(i) = 0.
      ra(:) = store(es_rx:es_rz,i)
      ha = store(es_h,i)
      ma = store(es_m,i)
      da = store(es_den,i)

      call getneighbours(i, nlistb, t0)
      tneib = tneib + t0
      do lb = 1, size(nlistb)
        j = nlistb(lb)
        rb(:) = store(es_rx:es_rz, j)
        mb = store(es_m, j)
        db = store(es_den, j)

        r(:) = ra(:)-rb(:)
        dr = sqrt(dot_product(r(:),r(:)))
        call get_w(dr, ha, w)
        dennew(i) = dennew(i) + mb * w
        currinterr = currinterr + mb/db * w
      end do
      call get_w(0., ha, w)
      dennew(i) = dennew(i) + ma * w
      store(es_h,i) = hfac * (ma / dennew(i))**(1./dim)
    end do
    !$omp end parallel do
    store(es_den,1:realpartnumb) = dennew(:)

    call system_clock(finish)
    call addTime(' circuit1', finish - start - tneib)
  end subroutine
end module
