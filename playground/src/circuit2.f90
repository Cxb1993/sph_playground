module circuit2
  use const
  use omp_lib
  use timing,           only: addTime
  use kernel,           only: get_hessian, &
                              get_krad, &
                              get_n2w, &
                              get_nw, &
                              get_hessian_rr
  use BC
  use neighboursearch,  only: getneighbours,&
                              getNeibListL1,&
                              getNeibListL2

  implicit none

  public :: c2init, c2

  private
  save
    integer(8)  :: start=0, finish=0
    integer     :: s_dim, s_ttp, s_ktp, s_dtp, s_adden, s_artts, initdone = 0
    real        :: s_kr
contains

  subroutine c2init()
    use state,  only: get_difftype, getdim, &
                      get_tasktype, get_kerntype, &
                      getAdvancedDensity, &
                      getArtificialTerms

    call getdim(s_dim)
    call get_krad(s_kr)
    call get_tasktype(s_ttp)
    call get_kerntype(s_ktp)
    call get_difftype(s_dtp)
    call getAdvancedDensity(s_adden)
    call getArtificialTerms(s_artts)
    initdone = 1
  end subroutine

  subroutine c2(c, ptype, pos, v, dv, mas, den, h, om, P, u, du, dh, cf, dcf, kcf, dfdx)
    real, allocatable, intent(in)    :: pos(:,:), v(:,:), mas(:), h(:), den(:), P(:), c(:),&
                                        u(:), cf(:,:), kcf(:,:,:), om(:)
    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(inout) :: dv(:,:), du(:), dh(:), dcf(:,:), dfdx(:,:,:)

    real                 :: dr, rhoa, rhob, qa(3), qb(3), qc, n2wa, r2, &
                            nwa(3), nwb(3), rab(3), vab(3), vba(3), urab(3), Pa(3), Pb(3), &
                            Hesa(3,3), Hesb(3,3), oddi, oddj, kcfij(3,3), ktmp
    integer, allocatable :: nlista(:), nlistb(:)
    integer              :: i, j, la, lb, n
    integer(8)           :: t0, tneib

    if (initdone == 0) then
      call c2init()
    end if
    call system_clock(start)
    n = size(ptype)
    tneib = 0.


    if (( s_ktp == 3 ).and.( s_dtp == 1 )) then
      call gradf(s_dim, mas, den, pos, v, h, om, dfdx)
    end if
    call getNeibListL1(nlista)
    !$omp parallel do default(none)&
    !$omp private(rab, dr, vab, urab, rhoa, rhob, nwa, nwb, qa, qb, qc, Pa, Pb)&
    !$omp private(n2wa, j, i, r2, oddi ,oddj, la, lb)&
    !$omp private(nlistb, Hesa, Hesb, vba, t0, kcfij, ktmp) &
    !$omp shared(dv, du, dh, dcf, n, pos, h, v, den, c, p, om, mas, u, kcf, cf)&
    !$omp shared(ptype, dfdx, nlista)&
    !$omp shared(s_dim, s_kr, s_ktp, s_dtp, s_ttp, s_adden, s_artts)&
    !$omp reduction(+:tneib)
    do la = 1, size(nlista)
      i = nlista(la)
      dv(:,i) = 0.
      du(i) = 0.
      dh(i) = 0.
      dcf(:,i) = 0.
      call getneighbours(i, pos, h, nlistb, t0)
      tneib = tneib + t0
      do lb = 1, size(nlistb)
        j = nlistb(lb)
        rab(:) = pos(:,i) - pos(:,j)
        r2 = dot_product(rab(:),rab(:))
        dr = sqrt(r2)
        vab(:) = v(:,i) - v(:,j)
        vba(:) = v(:,j) - v(:,i)
        urab(:) = rab(:) / dr
        select case (s_ttp)
        case (1, 9)
          ! hydroshock ! soundwave
          qa(:) = 0.
          qb(:) = 0.
          qc = 0.
          rhoa = den(i)
          rhob = den(j)
          call get_nw(rab, h(i), nwa)
          call get_nw(rab, h(j), nwb)
          if (s_artts == 1) then
            call art_viscosity(rhoa, rhob, vab, urab, rab, dr, c(i), c(j), om(i), om(j), h(i), h(j), qa, qb)
            call art_termcond(nwa, nwb, urab, P(i), P(j), u(i), u(j), rhoa, rhob, qc)
          end if
          ! Pa(:) = (P(i) + qa) * nwa(:) / (rhoa**2 * om(i))
          ! Pb(:) = (P(j) + qb) * nwb(:) / (rhob**2 * om(j))

          dv(:,i) = dv(:,i) - mas(j) * ( &
                      P(i) * nwa(:) / (rhoa**2 * om(i)) &
                      + P(j) * nwb(:) / (rhob**2 * om(j)) &
                      + qa(:) &
                      + qb(:) &
                    )

          du(i)   = du(i) + mas(j) * ( &
                      dot_product(vab(:), P(i) * nwa(:) / (rhoa**2 * om(i))) &
                      + dot_product(vab(:),qa(:)) &
                      + qc &
                    )

          if ( s_adden == 1 ) then
            dh(i) = dh(i) + mas(j) * dot_product(vab(:), nwa(:))
          end if
        case (2, 3)
          ! heatconduction
          ! call get_n2w(rab, h(i), n2wa)
          ! dcf(:,i) = dcf(:,i) - mas(j) / (den(i) * den(j)) * &
          !             2. * kcf(:,i) * kcf(:,j) / (kcf(:,i) + kcf(:,j)) * &
          !             (cf(:,i) - cf(:,j)) * n2wa

          ! kcfij(:,:) = 2. * kcf(:,:,i) * kcf(:,:,j) / (kcf(:,:,i) + kcf(:,:,j))
          call get_hessian(rab, h(i), Hesa)
          ! call get_hessian(-rab, h(j), Hesb)
          kcfij(:,:) = 0.
          ktmp = kcf(1,1,i) + kcf(1,1,j)
          if (abs(ktmp) > eps0) then
            kcfij(1,1) = kcf(1,1,i) * kcf(1,1,j) / ktmp
          end if
          ktmp = kcf(1,2,i) + kcf(1,2,j)
          if (abs(ktmp) > eps0) then
            kcfij(1,2) = kcf(1,2,i) * kcf(1,2,j) / ktmp
          end if
          ktmp = kcf(1,3,i) + kcf(1,3,j)
          if (abs(ktmp) > eps0) then
            kcfij(1,3) = kcf(1,3,i) * kcf(1,3,j) / ktmp
          end if
          kcfij(2,1) = kcfij(1,2)
          ktmp = kcf(2,2,i) + kcf(2,2,j)
          if (abs(ktmp) > eps0) then
            kcfij(2,2) = kcf(2,2,i) * kcf(2,2,j) / ktmp
          end if
          ktmp = kcf(2,3,i) + kcf(2,3,j)
          if (abs(ktmp) > eps0) then
            kcfij(2,3) = kcf(2,3,i) * kcf(2,3,j) / ktmp
          end if
          kcfij(3,1) = kcfij(1,3)
          kcfij(3,2) = kcfij(3,2)
          ktmp = kcf(3,3,i) + kcf(3,3,j)
          if (abs(ktmp) > eps0) then
            kcfij(3,3) = kcf(3,3,i) * kcf(3,3,j) / ktmp
          end if
          kcfij(:,:) = 2. * kcfij(:,:)

          ! Hesa(:,:) = 0.5*(Hesa(:,:)/den(i) + Hesb(:,:)/den(j))

          dcf(1,i) = dcf(1,i) + mas(j)/den(j) * (cf(1,j) - cf(1,i)) * &
            ( dot_product(kcfij(1,:),Hesa(1,:)) + &
              dot_product(kcfij(2,:),Hesa(2,:)) + &
              dot_product(kcfij(3,:),Hesa(3,:)) )
          ! if (abs(cf(1,i) - cf(1,j)) > epsilon(0.)) then
          !   print*, dr, h(i), mas(j)/(den(i) * den(j))
          !   print*,'--------------'
          !   print*, kcfij(1,:)
          !   print*, kcfij(2,:)
          !   print*, kcfij(3,:)
          !   print*,'--------------'
          !   print*, Hes(1,:)
          !   print*, Hes(2,:)
          !   print*, Hes(3,:)
          !   print*,'--------------'
          !   print*, 'dT_', i, ' = ',  dcf(1,i)
          !   print*,'--------------'
          !   read*
          ! end if
        case(4)
        case(5)
          ! 'diff-laplace'
          call get_n2w(rab, h(i), n2wa)
          dv(:,i)  = dv(:,i) + mas(j)/den(j) * vba(:) * n2wa
          ! call get_hessian(rab, h(i), Hes)
          ! kcfij(:) = 2. * kcf(:,i) * kcf(:,j) / (kcf(:,i) + kcf(:,j))
          ! dcf(1,i) = dcf(1,i) - mas(j) / (den(i) * den(j)) * (cf(1,i) - cf(1,j)) * &
          !           (kcfij(1)*Hes(1,1) + kcfij(2)*Hes(2,2) + kcfij(3)*Hes(3,3))
        case(6)
          ! diff-graddiv
          if ( s_dtp == 1) then
            ! diff form
            ! if (ktp /= 3) then
              ! n2w fab
              call get_hessian(rab, h(i), Hesa)
              dv(1,i) = dv(1,i) + mas(j)/den(j) * dot_product(vba,Hesa(:,1))
              dv(2,i) = dv(2,i) + mas(j)/den(j) * dot_product(vba,Hesa(:,2))
              dv(3,i) = dv(3,i) + mas(j)/den(j) * dot_product(vba,Hesa(:,3))
            ! else
            !   ! 2nw
            !   call get_nw(rab, h(i), nwa)
            !   ! qa = dfdx(1,1,i) + dfdx(2,2,i) + dfdx(3,3,i)
            !   ! qb = dfdx(1,1,j) + dfdx(2,2,j) + dfdx(3,3,j)
            !   qa = v(1,i)
            !   qb = v(1,j)
            !   dv(:,i) = dv(:,i) + mas(j) * (qb - qa) * nwa(:)
            ! end if
          elseif ( s_dtp == 2) then
            ! symm form
            if ( s_ktp /= 3 ) then
              ! n2w fab
            else
              ! 2nw
              call get_nw(rab, h(i), nwa)
              call get_nw(rab, h(j), nwb)
              oddi = 1./om(i)/den(i)/den(i)
              oddj = 1./om(j)/den(j)/den(j)
              ! qa = dfdx(1,1,i) + dfdx(2,2,i) + dfdx(3,3,i)
              ! qb = dfdx(1,1,j) + dfdx(2,2,j) + dfdx(3,3,j)
              qa = v(1,i)
              qb = v(1,j)
              dv(:,i) = dv(:,i) + mas(j) * (qa*nwa(:)*oddi + qb*nwb(:)*oddj)
            end if
          else
            print *, 'Diff type is not set in circuit2'
            stop
          end if
        case(10)
          ! diff-artvisc
          call get_hessian_rr(rab, h(i), Hesa)
          dv(1,i) = dv(1,i) - mas(j)/den(j) * dot_product(vab, Hesa(:,1))
          dv(2,i) = dv(2,i) - mas(j)/den(j) * dot_product(vab, Hesa(:,2))
          dv(3,i) = dv(3,i) - mas(j)/den(j) * dot_product(vab, Hesa(:,3))
        case default
          print *, 'Task type was not defined in circuit2.f90: line 240.'
          stop
        end select
      end do

      select case (s_ttp)
      case(1,2,3,4,9)
        if ( s_adden == 1 ) then
          dh(i) =  (- h(i) / (s_dim * den(i))) * dh(i) / om(i)
        end if
      case(5, 6, 10)
        ! diff-graddiv ! diff-laplace ! diff-artvisc
        if ( s_ktp == 3 ) then
          dv(:,i) = dv(:,i) * den(i)
        end if
      case default
        print *, 'Task type was not set in circuit2 outside circle'
        stop
      end select
      ! print*, i, den(i)
    end do
    !$omp end parallel do
    ! print*, du(1:5)
    call system_clock(finish)
    call addTime(' circuit2', finish - start - tneib)
  end subroutine

  pure subroutine art_termcond(nwa, nwb, urab, pa, pb, ua, ub, da, db, qc)
    real, intent(in)  :: pa, pb, da, db, ua, ub, &
                          nwa(3), nwb(3), urab(3)
    real, intent(out) :: qc
    real              :: vsigu

    vsigu = sqrt(abs(pa - pb)/(0.5*(da + db)))
    qc = vsigu * (ua - ub) / (0.5*(da + db)) * 0.5 * dot_product((nwa(:) + nwb(:)),urab(:))
  end subroutine

  pure subroutine art_viscosity(da, db, vab, &
      urab, rab, dr, &
      ca, cb, oa, ob, ha, hb, &
      qa, qb)
    real, intent(in)  :: da, db, vab(3), urab(3), rab(3),&
                         dr, ca, cb, oa, ob, ha, hb
    real, intent(out) :: qa(3), qb(3)
    real              :: alpha, betta, dvr, Hesrr(3,3), n2w(3)
    qa(:) = 0.
    qb(:) = 0.
    alpha = 1.
    betta = 2.

    ! qa * nwa(:) / (rhoa**2 * om(i))

    dvr = dot_product(vab,urab)
    if (dvr <= 0) then
      ! call get_n2w(rab, ha, n2w(1))
      call get_hessian_rr(rab, ha, Hesrr)
      qa(1) = -0.5 * (alpha*ca - betta*dvr) * &
          dot_product(vab, Hesrr(:,1)) * dr / (-2) / (da * oa)
      qa(2) = -0.5 * (alpha*ca - betta*dvr) * &
          dot_product(vab, Hesrr(:,2)) * dr / (-2) / (da * oa)
      qa(3) = -0.5 * (alpha*ca - betta*dvr) * &
          dot_product(vab, Hesrr(:,3)) * dr / (-2) / (da * oa)

      call get_hessian_rr(rab, hb, Hesrr)
      qb(1) = -0.5 * (alpha*cb - betta*dvr) * &
          dot_product(vab, Hesrr(:,1)) * dr / (-2) / (db * ob)
      qb(2) = -0.5 * (alpha*cb - betta*dvr) * &
          dot_product(vab, Hesrr(:,2)) * dr / (-2) / (db * ob)
      qb(3) = -0.5 * (alpha*cb - betta*dvr) * &
          dot_product(vab, Hesrr(:,3)) * dr / (-2) / (db * ob)
      ! call get_n2w(rab, hb, n2w(1))
      ! qb(:) = -0.5 * (alpha*cb - betta*dvr) * &
      !     vab(:) * n2w(1) * dr / (-2) / (db * ob)
    end if
  end subroutine

  subroutine gradf(dim, m, d, x, v, h, om, nv)
    real, allocatable, intent(in)  :: m(:), d(:), x(:,:), v(:,:), h(:), om(:)
    real, allocatable, intent(inout) :: nv(:,:,:)
    integer, intent(in)            :: dim

    integer, allocatable :: nlista(:), nlistb(:)
    integer              :: n, i, j, la, lb, ni, nj
    integer(8)           :: t0, tneib

    real    :: vba(3), rab(3), nwi(3), nwj(3), oddi, oddj

    n = size(m)
    ! if ( .not.allocated(nv) ) then
    !   allocate(nv(3,3,n))
    ! end if
    tneib = 0.
    ! call getNeibListL2(nlista)
    call getNeibListL1(nlista)

    !$omp parallel do default(none)&
    !$omp private(rab, vba, i, j, la, lb, ni, nj, nlistb, t0) &
    !$omp private(nwi, nwj, oddi, oddj)&
    !$omp shared(dim, m, d, x, v, h, om, nv, n, nlista)&
    !$omp reduction(+:tneib)
    do la = 1,size(nlista)
      i = nlista(la)
      nv(:,:,i) = 0.
      call getneighbours(i, x, h, nlistb, t0)
      tneib = tneib + t0
      do lb = 1, size(nlistb)
        j = nlistb(lb)
        vba(:) = v(:,j) - v(:,i)
        rab(:) = x(:,i) - x(:,j)
        call get_nw(rab, h(i), nwi)
        call get_nw(rab, h(j), nwj)
        oddi = 1./om(i)/d(i)/d(i)
        oddj = 1./om(j)/d(j)/d(j)
        do ni = 1,dim
          do nj = 1,dim
            ! nv(ni,nj,i) = nv(ni,nj,i) + m(j)/d(j)*vba(ni)*nwi(nj)
            nv(ni,nj,i) = nv(ni,nj,i) + m(j)*(v(ni,i)*nwi(nj)*oddi + v(ni,j)*nwj(nj)*oddj)
          end do
        end do
      end do
      nv(:,:,i) = nv(:,:,i) / om(i) / d(i)
    end do
    !$omp end parallel do
    call addTime(' circuit2', -tneib)
  end subroutine
end module
