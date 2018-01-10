module circuit1
  use const
  use omp_lib
  use timing,           only: addTime
  use state,            only: getdim, &
                              getddwtype, &
                              gcoordsys
  use kernel,           only: get_krad, &
                              get_dw_dh, &
                              nw, &
                              get_w
  use neighboursearch,  only: getneighbours,&
                              getNeibListL1,&
                              getNeibListL2

  implicit none

  public :: c1_init, c1, destroy

  private
  save
    real, allocatable :: slnint(:), resid(:), dennew(:)
    integer(8) :: start=0, finish=0

contains
  subroutine c1_init(n)
    integer, intent(in) :: n
    allocate(slnint(n))
    allocate(resid(n))
    allocate(dennew(n))
  end subroutine c1_init

  subroutine c1(ptype, pos, mas, hfac, h, den, om, cf, dcf, kcf)
    use state, only: getAdvancedDensity

    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(in)    :: pos(:,:), mas(:), cf(:,:)
    real, allocatable, intent(inout) :: h(:), den(:), om(:), dcf(:,:), kcf(:,:,:)
    real, intent(in)                 :: hfac

    integer :: ad

    call getAdvancedDensity(ad)

    if (ad == 1) then
      call c1advanced(ptype, pos, mas, hfac, h, den, om, cf, dcf, kcf)
    else
      call c1simple(ptype, pos, mas, hfac, h, den, cf, dcf, kcf)
      om(:) = 1.
    end if
  end subroutine

  subroutine destroy()
    deallocate(slnint)
    deallocate(resid)
    deallocate(dennew)
  end subroutine

  subroutine c1advanced(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
    use state,  only: ginitvar
    use BC,     only: getCrossRef, needcrosref

    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(in)    :: pos(:,:), mas(:), cf(:,:), kcf(:,:,:)
    real, allocatable, intent(inout) :: h(:), den(:), om(:), dcf(:,:)
    real, intent(in)     :: sk
    real                 :: w, dwdh, r(3), dr, r2, dfdh, fh, hn, vba(3), nwa(3), nwb(3)
    real                 :: allowerror, maxinterr, currinterr
    integer              :: n, i, rj, pj, la, lb, dim, iter, ktp
    integer(8)           :: t0, tneib
    integer, allocatable :: nlista(:), nlistb(:)

    call system_clock(start)

    ! n = size(den)

    call getdim(dim)
    call getddwtype(ktp)

    call getNeibListL1(nlista)
    allowerror = 1e-8
    slnint(:) = h(:)
    dennew(:) = den(:)
    resid(:)  = 0.
    resid(nlista) = 1.
    iter = 0
    tneib = 0.
    maxinterr  = 0.
    currinterr = 0.
    ! print*, nlista
    ! print*, '1'
    do while ((maxval(resid, mask=(resid>0)) > allowerror) .and. (iter < 100))
      maxinterr  = 0.
      iter = iter + 1
      !$omp parallel do default(none)&
      !$omp private(r, dr, dwdh, w, dfdh, fh, hn, rj, pj, i, la, lb, r2, t0, nlistb)&
      !$omp private(nwa, nwb, vba, currinterr)&
      !$omp shared(resid, allowerror, n, pos, mas, dim, sk, h, ktp, maxinterr)&
      !$omp shared(nlista, den, dennew, om, slnint, dcf, cf, needcrosref)&
      !$omp shared(nw)&
      !$omp reduction(+:tneib)
      do la = 1, size(nlista)
        i = nlista(la)
        if (resid(i) > allowerror) then
          dennew(i) = 0.
          om(i)  = 0.
          currinterr = 0.
          dcf(:,i) = 0.
          ! print*, i
          call getneighbours(i, pos, slnint, nlistb, t0)
          tneib = tneib + t0
          ! print*, nlista
          ! print*, "a=", i, "   neibs=", nlistb, pos(:,i)
          ! read*
          do lb = 1, size(nlistb)
            pj = nlistb(lb)                 ! phantom
            if (needcrosref == 1) then
              rj = getCrossRef(nlistb(lb))  ! real
            else
              rj = pj
            end if
            ! print*, "b=",lb, "id=",nlistb(lb), "real=", j
            ! read*
            r(:) = pos(:,i) - pos(:,pj)
            r2 = dot_product(r(:),r(:))
            ! print*, -3
            dr = sqrt(r2)
            ! print*, -2
            call get_dw_dh(dr, slnint(i), dwdh)
            ! print*, -1
            call get_w(dr, slnint(i), w)
            ! print*,  0
            dennew(i) = dennew(i) + mas(rj) * w
            om(i) = om(i) + mas(rj) * dwdh
            currinterr = currinterr + mas(rj)/den(rj) * w
            if (ktp == esd_2nw) then
              call nw(r(:), pos(:,i), pos(:,pj), dr, slnint(i), nwa)
              dcf(:,i) = dcf(:,i) + mas(rj)/den(rj)*(cf(1,rj) - cf(1,i))*nwa(:)
            end if
          end do
          ! ---------------------------------------------------------!
          !      There is no particle itself in neighbour list       !
          ! ---------------------------------------------------------!
          ! print*,'c1', 1
          call get_w(0., slnint(i), w)
          dennew(i) = dennew(i) + mas(i) * w
          currinterr = currinterr + mas(i)/den(i) * w
          if (currinterr > maxinterr) then
            maxinterr = currinterr
          end if

          call get_dw_dh(0., slnint(i), dwdh)
          om(i) = om(i) + mas(i) * dwdh
          ! -(**)----------------------------------------------------!
          ! print*,'c1', 4, om(i), mas(i), dwdh
          om(i) = 1. - om(i) * (-slnint(i) / (dim * dennew(i)))
          ! print*,'c1', 5, i, om(i), slnint(i), dim, den(i)
          dfdh = - dim * dennew(i) * om(i) / slnint(i)
          ! print*,'c1', 7, den(i), om(i), slnint(i)
          fh  = mas(i) * (sk / slnint(i)) ** dim - dennew(i)
          ! print*,'c1', 8, dfdh
          hn = slnint(i) - fh / dfdh
          if (hn <= 0.) then
            error stop "Negative smoothing length in c1 advanced"
          end if
          ! print*,'c1', 9
          resid(i) = abs(hn - slnint(i)) / h(i)
          slnint(i) = hn
        end if
        ! print *,'c1', dennew(i), slnint(i), om(i)
        ! read*
      end do
      !$omp end parallel do

      den(:) = dennew(:)
    end do
    ! print*, '2'
    if (iter > 10) then
      print*, "Warn: density NR: solution took ", iter, "iterations, with max norm error", maxval(resid, mask=(resid>0))
    end if

    ! if (abs(1. - maxinterr) > 0.1) then
    !   print*, "Warn: density NR: kernel integral condition does not fulfilled Int(V_{ab}*W_{ab}) = ", maxinterr
    ! end if
    ! read*
    h(:) = slnint(:)
    call system_clock(finish)
    call addTime(' circuit1', finish - start - tneib)
  end subroutine

! Direct density summation
  subroutine c1simple(ptype, pos, mas, sk, sln, den, cf, dcf, kcf)
    use BC, only: getCrossRef

    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(in)    :: pos(:,:), mas(:), cf(:,:), kcf(:,:,:)
    real,              intent(in)    :: sk
    real, allocatable, intent(inout) :: den(:), sln(:), dcf(:,:)
    real                             :: w, r(3), dr, currinterr
    integer                          :: i, rj, pj, la, lb, dim
    integer, allocatable             :: nlista(:), nlistb(:)
    integer(8)                       :: t0, tneib

    call system_clock(start)
    call getNeibListL2(nlista)
    call getdim(dim)

    dennew(:) = den(:)
    currinterr = 0.

    tneib = 0.
    !$omp parallel do default(none)&
    !$omp private(r, dr, w, rj, pj, i, la, lb, nlistb, t0, currinterr)&
    !$omp shared(pos, mas, sk, sln, dim)&
    !$omp shared(nlista, den, dennew)&
    !$omp reduction(+:tneib)
    do la = 1, size(nlista)
      i = nlista(la)

      dennew(i) = 0.
      ! currinterr = 0.
      call getneighbours(i, pos, sln, nlistb, t0)
      tneib = tneib + t0
      do lb = 1, size(nlistb)
        rj = getCrossRef(nlistb(lb))
        pj = nlistb(lb)
        r(:) = pos(:,i) - pos(:,pj)
        dr = sqrt(dot_product(r(:),r(:)))
        call get_w(dr, sln(i), w)
        dennew(i) = dennew(i) + mas(rj) * w
        currinterr = currinterr + mas(rj)/den(rj) * w
      end do
      call get_w(0., sln(i), w)
      dennew(i) = dennew(i) + mas(i) * w
      ! currinterr = currinterr + mas(i)/den(i) * w
      sln(i) = sk * (mas(i) / den(i))**(1./dim)
    end do
    !$omp end parallel do
    den(:) = dennew(:)

    call system_clock(finish)
    call addTime(' circuit1', finish - start - tneib)
  end subroutine
end module
