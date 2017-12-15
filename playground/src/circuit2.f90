module circuit2
  use const
  use omp_lib
  use timing,           only: addTime
  use kernel,           only: hessian, &
                              get_krad, &
                              n2w, &
                              nw, &
                              hessian_rr

  use BC
  use neighboursearch,  only: getneighbours,&
                              getNeibListL1,&
                              getNeibListL2

  implicit none

  public :: c2init, c2, c15

  private
  save
    integer(8)  :: start=0, finish=0
    integer     :: s_dim, s_ttp, s_ktp, s_adden, s_artts, s_ivt, initdone = 0, s_origin
    real        :: s_kr
    real, allocatable :: dcftmp(:,:)

contains

  subroutine c2init()
    use state,  only: get_difftype, getdim, &
                      get_tasktype, getkerntype, &
                      getAdvancedDensity, &
                      getArtificialTerms, &
                      getpartnum, &
                      ginitvar, &
                      gorigin
    integer :: n

    call getdim(s_dim)
    call get_krad(s_kr)
    call get_tasktype(s_ttp)
    call getkerntype(s_ktp)
    call getAdvancedDensity(s_adden)
    call getArtificialTerms(s_artts)
    call ginitvar(s_ivt)
    call gorigin(s_origin)

    if (s_ktp == 3) then
      call getpartnum(n)
      allocate(dcftmp(3,n))
    end if
    initdone = 1
  end subroutine

  subroutine c2(c, ptype, pos, v, dv, mas, den, h, om, P, u, du, dh, cf, dcf, kcf)
    use state, only: diff_conductivity, diff_isotropic, mhd_magneticconstant

    real, allocatable, intent(in)    :: pos(:,:), v(:,:), mas(:), h(:), den(:), P(:), c(:),&
                                        u(:), cf(:,:), om(:)
    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(inout) :: dv(:,:), du(:), dh(:), dcf(:,:), kcf(:,:,:)

    real                 :: dr, rhoa, rhob, qa(3), qb(3), qc, qd(3), n2wa, r2, &
                            nwa(3), nwb(3), rab(3), vab(3), vba(3), urab(3), &
                            Hesa(3,3), odda, oddb, kcfa(3,3), kcfb(3,3), kcfab(3,3), ktmp, phi, &
                            tmpt1(3,3), tmpt2(3,3), tmpt3(3,3), Ma(3,3), Mb(3,3), Mc(3)

    integer, allocatable :: nlista(:), nlistb(:)
    integer              :: i, j, la, lb, n, li, lj
    integer(8)           :: t0, tneib

    call system_clock(start)

    if (initdone == 0) then
      call c2init()
    end if

    if (s_ktp == 3) then
      dcftmp(:,:) = dcf(:,:)
      dcf(:,:) = 0.
    end if

    n = size(ptype)
    tneib = 0.
    ! print*, 1
    call getNeibListL1(nlista)

    !$omp parallel do default(none)&
    !$omp private(rab, dr, vab, urab, rhoa, rhob, nwa, nwb, qa, qb, qc, qd)&
    !$omp private(n2wa, j, i, r2, odda ,oddb, la, lb)&
    !$omp private(nlistb, Hesa, vba, t0, kcfa, kcfb, kcfab, ktmp, phi) &
    !$omp private(tmpt1, tmpt2, tmpt3, li, lj, Ma, Mb, Mc)&
    !$omp shared(dv, du, dh, dcf, n, pos, h, v, den, c, p, om, mas, u, kcf, cf)&
    !$omp shared(ptype, nlista, dcftmp)&
    !$omp shared(s_dim, s_kr, s_ktp, s_ttp, s_adden, s_artts, s_ivt)&
    !$omp shared(nw, hessian, hessian_rr)&
    !$omp shared(diff_conductivity, diff_isotropic, mhd_magneticconstant)&
    !$omp reduction(+:tneib)
    overa: do la = 1, size(nlista)
      i = nlista(la)
      dv(:,i) = 0.
      du(i) = 0.
      dh(i) = 0.
      dcf(:,i) = 0.
      kcf(:,2,i) = 0.
      tmpt1(:,:) = 0.
      tmpt2(:,:) = 0.
      tmpt3(:,:) = 0.
      call getneighbours(i, pos, h, nlistb, t0)
      tneib = tneib + t0
      select case (s_ttp)
      case (eeq_diffusion, eeq_magnetohydrodiffusion)
        kcfa(:,:) = 0.
        kcfb(:,:) = 0.
        if (diff_isotropic > 0.) then
          kcfa(1,1) = diff_conductivity
          kcfa(2,2) = diff_conductivity
          kcfa(3,3) = diff_conductivity
        else
          do li = 1, 3
            do lj = 1,3
              kcfa(li,lj) = diff_conductivity*kcf(li,1,i)*kcf(lj,1,i)
            end do
          end do
        end if
      end select
      overb: do lb = 1, size(nlistb)
        j = nlistb(lb)
        qa(:) = 0.
        qb(:) = 0.
        qc    = 0.
        qd(:) = 0.
        rab(:) = pos(:,i) - pos(:,j)
        r2 = dot_product(rab(:),rab(:))
        dr = sqrt(r2)
        vab(:) = v(:,i) - v(:,j)
        vba(:) = v(:,j) - v(:,i)
        urab(:) = rab(:) / dr
        select case (s_ttp)
        case (eeq_hydro)
          rhoa = den(i)
          rhob = den(j)
          call nw(rab, pos(:,i), pos(:,j), dr, h(i), nwa)
          call nw(rab, pos(:,i), pos(:,j), dr, h(j), nwb)

          if (s_artts == 1) then
            call art_viscosity(rhoa, rhob, vab, urab, rab, dr, s_dim, c(i), c(j), om(i), om(j), h(i), h(j), qa, qb)
            call art_termcond(nwa, nwb, urab, P(i), P(j), u(i), u(j), rhoa, rhob, om(i), om(j), qc)
          end if

          dv(:,i) = dv(:,i) - mas(j) * ( &
                      (P(i) * nwa(:) + qa(:)) / (rhoa**2 * om(i)) + &
                      (P(j) * nwb(:) + qb(:)) / (rhob**2 * om(j)) &
                    )

          du(i)   = du(i) + mas(j) * ( &
                      dot_product(vab(:), P(i) * nwa(:)) / (rhoa**2 * om(i)) - &
                      dot_product(vab(:), qa(:))/(rhoa * om(i)) + &
                      qc &
                    )

          if ( s_adden == 1 ) then
            dh(i) = dh(i) + mas(j) * dot_product(vab(:), nwa(:))
          end if
        case (eeq_magnetohydro)
          Ma(:,:) = 0.
          Mb(:,:) = 0.
          Mc(:)   = 0.
          rhoa = den(i)
          rhob = den(j)
          call nw(rab, pos(:,i), pos(:,j), dr, h(i), nwa)
          call nw(rab, pos(:,i), pos(:,j), dr, h(j), nwb)

          if (s_artts == 1) then
            call art_viscosity(rhoa, rhob, vab, urab, rab, dr, &
                                s_dim, c(i), c(j), om(i), om(j), h(i), h(j), qa, qb)
            call art_termcond(nwa, nwb, urab, P(i), P(j), u(i), u(j), rhoa, rhob, om(i), om(j), qc)
            call art_fdivbab(kcf(:,1,i), kcf(:,1,j), mas(j), nwa, nwb, om(i), om(j), den(i), den(j), qd)
          end if

          Mc(1) = dot_product(kcf(:,1,i),kcf(:,1,i))
          Mc(2) = dot_product(kcf(:,1,j),kcf(:,1,j))
          do li = 1,3
            do lj = 1,3
              if (li == lj) then
                Ma(li,lj) = P(i) + (0.5*Mc(1) - kcf(li,1,i)*kcf(lj,1,i)) /mhd_magneticconstant
                Mb(li,lj) = P(j) + (0.5*Mc(2) - kcf(li,1,j)*kcf(lj,1,j)) /mhd_magneticconstant
              else
                Ma(li,lj) = -kcf(li,1,i)*kcf(lj,1,i) /mhd_magneticconstant
                Mb(li,lj) = -kcf(li,1,j)*kcf(lj,1,j) /mhd_magneticconstant
              end if
            end do
          end do

          Ma(1,1) = dot_product(Ma(1,:), nwa(:)) + qa(1)
          Ma(2,1) = dot_product(Ma(2,:), nwa(:)) + qa(2)
          Ma(3,1) = dot_product(Ma(3,:), nwa(:)) + qa(3)

          Mb(1,1) = dot_product(Mb(1,:), nwb(:)) + qb(1)
          Mb(2,1) = dot_product(Mb(2,:), nwb(:)) + qb(2)
          Mb(3,1) = dot_product(Mb(3,:), nwb(:)) + qb(3)

          dv(:,i) = dv(:,i) - mas(j) * ( &
                      Ma(:,1) / (rhoa**2 * om(i)) + &
                      Mb(:,1) / (rhob**2 * om(j)) &
                    ) - qd(:)

          du(i)   = du(i) + mas(j) * ( &
                      dot_product(vab(:), nwa(:))*P(i)/(rhoa**2 * om(i)) - &
                      dot_product(vab(:), qa(:))/(rhoa * om(i)) + &
                      qc &
                    )

          kcf(:,2,i) = kcf(:,2,i) - 1./(rhoa * om(i))*mas(j)*( &
                      vab(:) * dot_product(kcf(:,1,i),nwa(:)) - &
                      kcf(:,1,i) * dot_product(vab(:),nwa(:)) &
                    )

          if ( s_adden == 1 ) then
            dh(i) = dh(i) + mas(j) * dot_product(vab(:), nwa(:))
          end if
        case (eeq_diffusion)
          if (diff_isotropic > 0.) then
            kcfb(1,1) = diff_conductivity
            kcfb(2,2) = diff_conductivity
            kcfb(3,3) = diff_conductivity
          else
            do li = 1, 3
              do lj = 1,3
                kcfb(li,lj) = diff_conductivity*kcf(li,1,j)*kcf(lj,1,j)
              end do
            end do
          end if

          call nw(rab, pos(:,i), pos(:,j), dr, h(i), nwa)
          call nw(rab, pos(:,i), pos(:,j), dr, h(j), nwb)

          if (s_ktp == 1) then
            call hessian(rab, pos(:,i), pos(:,j), h(i), Hesa)
            do li = 1, 3
              do lj = 1,3
                tmpt1(li,lj) = tmpt1(li,lj) + mas(j)/den(j) * &
                                (kcfb(li,lj) - kcfa(li,lj)) * nwa(li)
                tmpt2(li,lj) = tmpt2(li,lj) + mas(j)/den(j) * &
                                (cf(1,j) - cf(1,i)) * nwa(lj)
                tmpt3(li,lj) = tmpt3(li,lj) + kcfa(li,lj) *&
                                mas(j)/den(j) * (cf(1,j) - cf(1,i)) * Hesa(li,lj)
              end do
            end do
          else if ((s_ktp == 2).or.(s_ktp==4)) then
            kcfab(:,:) = (kcfa(:,:)+kcfb(:,:))/2.
            call hessian(rab, pos(:,i), pos(:,j), h(i), Hesa)
            dcf(1,i) = dcf(1,i) + mas(j)/den(j) * (cf(1,j) - cf(1,i)) * &
              ( dot_product(kcfab(1,:),Hesa(1,:)) + &
                dot_product(kcfab(2,:),Hesa(2,:)) + &
                dot_product(kcfab(3,:),Hesa(3,:)) )
          else if (s_ktp == 3) then
            odda = 1./om(i)/den(i)/den(i)
            oddb = 1./om(j)/den(j)/den(j)
            qa(1) = dot_product(kcfa(1,:),dcftmp(:,i))
            qa(2) = dot_product(kcfa(2,:),dcftmp(:,i))
            qa(3) = dot_product(kcfa(3,:),dcftmp(:,i))
            qb(1) = dot_product(kcfb(1,:),dcftmp(:,j))
            qb(2) = dot_product(kcfb(2,:),dcftmp(:,j))
            qb(3) = dot_product(kcfb(3,:),dcftmp(:,j))

            dcf(1,i) = dcf(1,i) + mas(j)*( &
              dot_product(qa(:),nwa(:))*odda + dot_product(qb(:),nwb(:))*oddb)
          end if
        case (eeq_magnetohydrodiffusion)
          Ma(:,:) = 0.
          Mb(:,:) = 0.
          Mc(:)   = 0.
          rhoa = den(i)
          rhob = den(j)
          call nw(rab, pos(:,i), pos(:,j), dr, h(i), nwa)
          call nw(rab, pos(:,i), pos(:,j), dr, h(j), nwb)
          if (diff_isotropic > 0.) then
            kcfb(1,1) = 1.
            kcfb(:,1) = [1.,0.,0.]
            kcfb(:,2) = [0.,1.,0.]
            kcfb(:,3) = [0.,0.,1.]
          else
            do li = 1, 3
              do lj = 1,3
                kcfb(li,lj) = kcf(li,1,j)*kcf(lj,1,j)
              end do
            end do
          end if
          kcfb(:,:) = diff_conductivity*kcfb(:,:)
          ! print*, 1
          ! print*, mhd_magneticconstant,
          if (s_artts == 1) then
            call art_viscosity(rhoa, rhob, vab, urab, rab, dr, &
                                s_dim, c(i), c(j), om(i), om(j), h(i), h(j), qa, qb)
            ! call art_termcond(nwa, nwb, urab, P(i), P(j), u(i), u(j), rhoa, rhob, om(i), om(j), qc)
            call art_fdivbab(kcf(:,1,i), kcf(:,1,j), mas(j), nwa, nwb, om(i), om(j), den(i), den(j), qd)
          end if
          ! print*, 2
          Mc(1) = dot_product(kcf(:,1,i),kcf(:,1,i))
          Mc(2) = dot_product(kcf(:,1,j),kcf(:,1,j))
          do li = 1,3
            do lj = 1,3
              if (li == lj) then
                Ma(li,lj) = P(i) + (0.5*Mc(1) - kcf(li,1,i)*kcf(lj,1,i)) /mhd_magneticconstant
                Mb(li,lj) = P(j) + (0.5*Mc(2) - kcf(li,1,j)*kcf(lj,1,j)) /mhd_magneticconstant
              else
                Ma(li,lj) = -kcf(li,1,i)*kcf(lj,1,i) /mhd_magneticconstant
                Mb(li,lj) = -kcf(li,1,j)*kcf(lj,1,j) /mhd_magneticconstant
              end if
            end do
          end do
          ! print*, 3
          ! print*, Ma(:,1)
          ! print*, Ma(:,2)
          ! print*, Ma(:,3)
          ! print*, '-------'
          ! print*, i, kcf(:,1,i)
          ! print*, j, kcf(:,1,j)
          ! print*,'---------'
          ! print*, i, j, "|", om(i), om(j)
          ! read*
          Ma(1,1) = dot_product(Ma(1,:), nwa(:)) + qa(1)
          Ma(2,1) = dot_product(Ma(2,:), nwa(:)) + qa(2)
          Ma(3,1) = dot_product(Ma(3,:), nwa(:)) + qa(3)

          Mb(1,1) = dot_product(Mb(1,:), nwb(:)) + qb(1)
          Mb(2,1) = dot_product(Mb(2,:), nwb(:)) + qb(2)
          Mb(3,1) = dot_product(Mb(3,:), nwb(:)) + qb(3)

          dv(:,i) = dv(:,i) - mas(j) * ( &
                      Ma(:,1) / (rhoa**2 * om(i)) + &
                      Mb(:,1) / (rhob**2 * om(j)) &
                    ) - qd(:)

          du(i)   = du(i) + mas(j) * ( &
                      dot_product(vab(:), nwa(:))*P(i)/(rhoa**2 * om(i)) - &
                      dot_product(vab(:), qa(:))/(rhoa * om(i)) + &
                      qc &
                    )

          kcf(:,2,i) = kcf(:,2,i) - 1./(rhoa * om(i))*mas(j)*( &
                      vab(:) * dot_product(kcf(:,1,i),nwa(:)) - &
                      kcf(:,1,i) * dot_product(vab(:),nwa(:)) &
                    )

          if ( s_adden == 1 ) then
            dh(i) = dh(i) + mas(j) * dot_product(vab(:), nwa(:))
          end if
          ! diffusion
          if (s_ktp == 1) then
            call hessian(rab, pos(:,i), pos(:,j), h(i), Hesa)
            do li = 1, 3
              do lj = 1,3
                tmpt1(li,lj) = tmpt1(li,lj) + mas(j)/den(j) * &
                                (kcfb(li,lj) - kcfa(li,lj)) * nwa(li)
                tmpt2(li,lj) = tmpt2(li,lj) + mas(j)/den(j) * &
                                (cf(1,j) - cf(1,i)) * nwa(lj)
                tmpt3(li,lj) = tmpt3(li,lj) + kcfa(li,lj) * &
                                mas(j)/den(j) * (cf(1,j) - cf(1,i)) * Hesa(li,lj)
              end do
            end do
          else if ((s_ktp == 2).or.(s_ktp==4)) then
            kcfab(:,:) = (kcfa(:,:)+kcfb(:,:))/2.
            call hessian(rab, pos(:,i), pos(:,j), h(i), Hesa)
            dcf(1,i) = dcf(1,i) + mas(j)/den(j) * (cf(1,j) - cf(1,i)) * &
              ( dot_product(kcfab(1,:),Hesa(1,:)) + &
                dot_product(kcfab(2,:),Hesa(2,:)) + &
                dot_product(kcfab(3,:),Hesa(3,:)) )
          else if (s_ktp == 3) then
            odda = 1./om(i)/den(i)/den(i)
            oddb = 1./om(j)/den(j)/den(j)
            qa(1) = dot_product(kcfa(1,:),dcftmp(:,i))
            qa(2) = dot_product(kcfa(2,:),dcftmp(:,i))
            qa(3) = dot_product(kcfa(3,:),dcftmp(:,i))
            qb(1) = dot_product(kcfb(1,:),dcftmp(:,j))
            qb(2) = dot_product(kcfb(2,:),dcftmp(:,j))
            qb(3) = dot_product(kcfb(3,:),dcftmp(:,j))

            dcf(1,i) = dcf(1,i) + mas(j)*( &
              dot_product(qa(:),nwa(:))*odda + dot_product(qb(:),nwb(:))*oddb)
          end if
        ! case(5)
        !   ! 'diff-laplace'
        !   call n2w(rab, h(i), n2wa)
        !   dv(:,i)  = dv(:,i) + mas(j)/den(j) * vba(:) * n2wa
        !   ! call get_hessian(rab, h(i), Hes)
        !   ! kcfij(:) = 2. * kcf(:,i) * kcf(:,j) / (kcf(:,i) + kcf(:,j))
        !   ! dcf(1,i) = dcf(1,i) - mas(j) / (den(i) * den(j)) * (cf(1,i) - cf(1,j)) * &
        !   !           (kcfij(1)*Hes(1,1) + kcfij(2)*Hes(2,2) + kcfij(3)*Hes(3,3))
        ! case(6)
        !   call hessian(rab, pos(:,i), pos(:,j), h(i), Hesa)
        !   dv(1,i) = dv(1,i) + mas(j)/den(j) * dot_product(vba,Hesa(:,1))
        !   dv(2,i) = dv(2,i) + mas(j)/den(j) * dot_product(vba,Hesa(:,2))
        !   dv(3,i) = dv(3,i) + mas(j)/den(j) * dot_product(vba,Hesa(:,3))
        ! case(10)
        !   ! diff-artvisc
        !   call hessian_rr(rab, h(i), Hesa)
        !   dv(1,i) = dv(1,i) - mas(j)/den(j) * dot_product(vab, Hesa(:,1))
        !   dv(2,i) = dv(2,i) - mas(j)/den(j) * dot_product(vab, Hesa(:,2))
        !   dv(3,i) = dv(3,i) - mas(j)/den(j) * dot_product(vab, Hesa(:,3))
        case default
          print *, 'Task type was not defined in circuit2.f90: line 240.'
          stop
        end select
      end do overb

      select case (s_ttp)
      case(eeq_hydro, eeq_magnetohydro)
        if ( s_adden == 1 ) then
          dh(i) =  (- h(i) / (s_dim * den(i))) * dh(i) / om(i)
        end if
      case(eeq_diffusion, eeq_magnetohydrodiffusion)
        if ( s_adden == 1 ) then
          dh(i) =  (- h(i) / (s_dim * den(i))) * dh(i) / om(i)
        end if
        if (s_ktp == 1) then
          do li = 1,3
            do lj = 1,3
              dcf(1,i) = dcf(1,i) + tmpt1(li,lj) * tmpt2(li,lj) + tmpt3(li,lj)
            end do
          end do
        end if
      ! case(5, 6, 10)
      !   ! diff-graddiv ! diff-laplace ! diff-artvisc
      !   if ( s_ktp == 3 ) then
      !     dv(:,i) = dv(:,i) * den(i)
      !   end if
      case default
        print *, 'Task type was not set in circuit2 outside circle'
        stop
      end select
      ! print*, i, den(i)
    end do overa
    !$omp end parallel do
    ! print*, 2
    call system_clock(finish)
    call addTime(' circuit2', finish - start - tneib)
  end subroutine

  pure subroutine art_termcond(nwa, nwb, urab, pa, pb, ua, ub, da, db, oa, ob, qc)
    real, intent(in)  :: pa, pb, da, db, ua, ub, oa, ob, &
                          nwa(3), nwb(3), urab(3)
    real, intent(out) :: qc
    real              :: vsigu
    qc = 0.
    vsigu = sqrt(abs(pa - pb)/(0.5*(da + db)))
    qc = vsigu * (ua - ub) * 0.5 * (dot_product((nwa(:)/oa/da + nwb(:)/ob/db),urab(:)))
  end subroutine

  pure subroutine art_viscosity(da, db, vab, &
      urab, rab, dr, dim, &
      ca, cb, oa, ob, ha, hb, &
      qa, qb)
    real, intent(in)  :: da, db, vab(3), urab(3), rab(3),&
                         dr, ca, cb, oa, ob, ha, hb
    integer,intent(in):: dim
    real, intent(out) :: qa(3), qb(3)
    real              :: alpha, beta, dvr, Hesrr(3,3), vsig

    qa(:) = 0.
    qb(:) = 0.
    alpha = 1.
    beta  = 2.

    dvr = dot_product(vab,urab)
    if (dvr < 0) then
      vsig = alpha*ca + beta*abs(dvr)
      call hessian_rr(rab, ha, Hesrr)
      qa(1) = 0.5 * da * vsig * dot_product(vab, Hesrr(:,1)) * dr / (dim + 2.)! * dvr
      qa(2) = 0.5 * da * vsig * dot_product(vab, Hesrr(:,2)) * dr / (dim + 2.)! * dvr
      qa(3) = 0.5 * da * vsig * dot_product(vab, Hesrr(:,3)) * dr / (dim + 2.)! * dvr

      vsig = alpha*cb + beta*abs(dvr)
      call hessian_rr(rab, hb, Hesrr)
      qb(1) = 0.5 * db * vsig * dot_product(vab, Hesrr(:,1)) * dr / (dim + 2.)! * dvr
      qb(2) = 0.5 * db * vsig * dot_product(vab, Hesrr(:,2)) * dr / (dim + 2.)! * dvr
      qb(3) = 0.5 * db * vsig * dot_product(vab, Hesrr(:,3)) * dr / (dim + 2.)! * dvr
    end if
  end subroutine

  ! subroutine art_fdivbab(Ba, Bb, mb, nwa, nwb, oa, ob, da, db, qd)
  pure subroutine art_fdivbab(Ba, Bb, mb, nwa, nwb, oa, ob, da, db, qd)
    real, intent(in)  :: Ba(3), Bb(3), mb, nwa(3), nwb(3), oa, ob, da, db
    real, intent(out) :: qd(3)

    qd(:) = Ba(:)*mb*(dot_product(Ba(:),nwa(:))/oa/da/da + dot_product(Bb(:),nwb(:))/ob/db/db)
  end subroutine art_fdivbab

  ! pure subroutine art_resistivity(nwa, nwb, urab, pa, pb, ua, ub, da, db, oa, ob, qbt)
  !   real, intent(in)  :: pa, pb, da, db, ua, ub, oa, ob, &
  !                         nwa(3), nwb(3), urab(3)
  !   real, intent(out) :: qbt
  !   real              :: vsigu
  !   qc = 0.
  !
  !   vsigB = sqrt((dvx - projv*runix)**2 + (dvy - projv*runiy)**2 + (dvz - projv*runiz)**2)
  !   dBdissterm = (avBterm*grkerni + avBtermj*grkernj)*vsigB
  ! end subroutine

  subroutine c15(pos, mas, h, den, cf, om, dcf)
    real, allocatable, intent(in)    :: pos(:,:), mas(:), h(:), den(:),&
                                        cf(:,:), om(:)
    real, allocatable, intent(inout) :: dcf(:,:)

    real                 :: nwa(3), nwb(3), rab(3), dr
    integer, allocatable :: nlista(:), nlistb(:)
    integer              :: i, j, la, lb
    integer(8)           :: t0, tneib

    call system_clock(start)
    tneib = 0

    call getNeibListL1(nlista)
    do la = 1, size(nlista)
      i = nlista(la)
      call getneighbours(i, pos, h, nlistb, t0)
      tneib = tneib + t0
      do lb = 1, size(nlistb)
        j = nlistb(lb)
        rab(:) = pos(:,i) - pos(:,j)
        dr = sqrt(dot_product(rab(:),rab(:)))
        call nw(rab, pos(:,i), pos(:,j), dr, h(i), nwa)
        call nw(rab, pos(:,i), pos(:,j), dr, h(j), nwb)
        dcf(:,i) = dcf(:,i) + mas(j)*( &
            cf(1,i)*nwa(:)/om(i)/den(i)/den(i) + &
            cf(1,j)*nwb(:)/om(j)/den(j)/den(j) )
      end do
    end do

    call system_clock(finish)
    call addTime(' circuit-1.5', finish - start - tneib)
  end subroutine
end module
