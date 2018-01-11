module circuit2
  use const
  use omp_lib
  use timing,           only: addTime
  use kernel,           only: hessian, &
                              get_krad, &
                              n2w, &
                              nw, &
                              hessian_rr
  use neighboursearch,  only: getneighbours,&
                              getNeibListL1,&
                              getNeibListL2
  use state,            only: get_difftype, getdim, &
                              get_tasktype, getddwtype, &
                              getAdvancedDensity, &
                              getArtificialTerms, &
                              getpartnum, &
                              ginitvar, &
                              gorigin, &
                              getdiffisotropic, &
                              getdiffconductivity, &
                              getmhdmagneticpressure

  implicit none

  public :: c2init, c2

  private
  save
    integer(8)  :: start=0, finish=0
    integer     :: s_dim, s_ttp, s_ktp, s_adden, s_artts, s_ivt, initdone = 0, s_origin
    real        :: s_kr
    real, allocatable :: dcftmp(:,:)

contains

  subroutine c2init()
    integer :: n

    call getdim(s_dim)
    call get_krad(s_kr)
    call get_tasktype(s_ttp)
    call getddwtype(s_ktp)
    call getAdvancedDensity(s_adden)
    call getArtificialTerms(s_artts)
    call ginitvar(s_ivt)
    call gorigin(s_origin)

    if (s_ktp == esd_2nw) then
      call getpartnum(n)
      allocate(dcftmp(3,n))
    end if
    initdone = 1
  end subroutine

  subroutine c2(store)
    use BC,     only: getCrossRef, needcrosref
    real, allocatable, intent(inout) :: store(:,:)

    real :: &
      dr, rhoa, rhob, ra(3), rb(3), &
      qa(3), qb(3), qc, qd(3), r2, &
      nwa(3), nwb(3), rab(3), vab(3), vba(3), urab(3), &
      Hesa(3,3), odda, oddb, kta(3,3), ktb(3,3), ktab(3,3), &
      tmpt1(3,3), tmpt2(3,3), tmpt3(3,3), MPa(3,3), MPb(3,3), Mc(3),&
      difcond, mhdmuzero,&
      dva(3), dua, dha, ddta, dba(3), ba(3), bb(3), va(3), vb(3),&
      ha, hb, ca, cb, pa, pb, ua, ub, oma, omb, ta, tb, ma, mb

    integer, allocatable :: nlista(:), nlistb(:)
    integer :: &
      i, rj, pj, la, lb, li, lj, difiso
    integer(8)           :: t0, tneib

    call system_clock(start)
    call getmhdmagneticpressure(mhdmuzero)
    call getdiffisotropic(difiso)
    call getdiffconductivity(difcond)

    if (initdone == 0) then
      call c2init()
    end if

    if (s_ktp == esd_2nw) then
      dcftmp(:,:) = store(es_dtdx:es_dtdz,:)
      store(es_dtdx:es_dtdz,:) = 0.
    end if

    tneib = 0.
    call getNeibListL1(nlista)
    !$omp parallel do default(none)&
    !$omp private(rab, dr, ra, rb, vab, urab, rhoa, rhob, nwa, nwb, qa, qb, qc, qd)&
    !$omp private(rj, pj, i, r2, odda ,oddb, la, lb)&
    !$omp private(nlistb, Hesa, vba, t0, kta, ktb, ktab)&
    !$omp private(tmpt1, tmpt2, tmpt3, li, lj, MPa, MPb, Mc)&
    !$omp private(ha, hb, ca, cb, pa, pb, ua, ub, oma, omb, ta, tb, ma, mb)&
    !$omp private(dva, dua, dha, ddta, dba, ba, bb, va, vb)&
    !$omp shared(store, nlista, dcftmp)&
    !$omp shared(s_dim, s_kr, s_ktp, s_ttp, s_adden, s_artts, s_ivt)&
    !$omp shared(nw, hessian, hessian_rr)&
    !$omp shared(difcond, difiso, mhdmuzero)&
    !$omp shared(needcrosref)&
    !$omp reduction(+:tneib)
    overa: do la = 1, size(nlista)
      i = nlista(la)

      rhoa  = store(es_den,i)
      ra(:) = store(es_rx:es_rz,i)
      ba(:) = store(es_bx:es_bz,i)
      va(:) = store(es_vx:es_vz,i)
      ha  = store(es_h,i)
      ca  = store(es_c,i)
      pa  = store(es_p,i)
      ua  = store(es_u,i)
      oma = store(es_om,i)
      ta  = store(es_t,i)
      ma  = store(es_m,i)

      dva(:) = 0.
      dua = 0.
      dha = 0.
      ddta = 0.
      dba(:) = 0.

      tmpt1(:,:) = 0.
      tmpt2(:,:) = 0.
      tmpt3(:,:) = 0.
      kta(:,:) = 0.
      ktb(:,:) = 0.

      call getneighbours(i, nlistb, t0)
      tneib = tneib + t0
      select case (s_ttp)
      case (eeq_diffusion, eeq_magnetohydrodiffusion)
        if (difiso == 1) then
          kta(:,1) = [difcond,0.,0.]
          kta(:,2) = [0.,difcond,0.]
          kta(:,3) = [0.,0.,difcond]

          ktb(:,1) = [difcond,0.,0.]
          ktb(:,2) = [0.,difcond,0.]
          ktb(:,3) = [0.,0.,difcond]
        else
          do li = 1, 3
            do lj = 1,3
              kta(li,lj) = difcond*ba(li)*ba(lj)
            end do
          end do
        end if
      end select
      overb: do lb = 1, size(nlistb)
        pj = nlistb(lb)
        if (needcrosref == 1) then
          rj = getCrossRef(nlistb(lb))  ! real
        else
          rj = pj
        end if

        rhob = store(es_den,rj)
        rb(:) = store(es_rx:es_rz,pj)
        bb(:) = store(es_bx:es_bz,rj)
        vb(:) = store(es_vx:es_vz,rj)
        hb = store(es_h,rj)
        cb = store(es_c,rj)
        pb = store(es_p,rj)
        ub = store(es_u,rj)
        omb = store(es_om,rj)
        tb = store(es_t,rj)
        mb = store(es_m,rj)

        qa(:) = 0.
        qb(:) = 0.
        qc    = 0.
        qd(:) = 0.
        rab(:) = ra(:) - rb(:)
        r2 = dot_product(rab(:),rab(:))
        dr = sqrt(r2)
        vab(:) = va(:) - vb(:)
        vba(:) = vb(:) - va(:)
        urab(:) = rab(:) / dr
        select case (s_ttp)
        case (eeq_hydro)
          call nw(rab, ra(:), rb(:), dr, ha, nwa)
          call nw(rab, ra(:), rb(:), dr, hb, nwb)

          if (s_artts == 1) then
            call art_viscosity(rhoa, rhob, vab, urab, rab, dr, s_dim, ca, cb, ha, hb, qa, qb)
            call art_termcond(nwa, nwb, urab, pa, pb, ua, ub, rhoa, rhob, oma, omb, qc)
          end if

          dva(:) = dva(:) - mb * ( &
                      (pa * nwa(:) + qa(:)) / (rhoa**2 * oma) + &
                      (pb * nwb(:) + qb(:)) / (rhob**2 * omb) &
                    )

          dua   = dua + mb * ( &
                      dot_product(vab(:), pa * nwa(:)) / (rhoa**2 * oma) - &
                      dot_product(vab(:), qa(:))/(rhoa * oma) + &
                      qc &
                    )

        case (eeq_magnetohydro)
          MPa(:,:) = 0.
          MPb(:,:) = 0.
          Mc(:)   = 0.
          call nw(rab, ra, rb, dr, ha, nwa)
          call nw(rab, ra, rb, dr, hb, nwb)

          if (s_artts == 1) then
            call art_viscosity(rhoa, rhob, vab, urab, rab, dr, &
                                s_dim, ca, cb, ha, hb, qa, qb)
            call art_termcond(nwa, nwb, urab, pa, pb, ua, ub, rhoa, rhob, oma, omb, qc)
            call art_fdivbab(ba, bb, ma, nwa, nwb, oma, omb, rhoa, rhob, qd)
          end if

          Mc(1) = dot_product(ba(:),ba(:))
          Mc(2) = dot_product(bb(:),bb(:))
          do li = 1,3
            do lj = 1,3
              if (li == lj) then
                MPa(li,lj) = pa + (0.5*Mc(1) - ba(li)*ba(lj)) /mhdmuzero
                MPb(li,lj) = pb + (0.5*Mc(2) - bb(li)*bb(lj)) /mhdmuzero
              else
                MPa(li,lj) = -ba(li)*ba(lj) /mhdmuzero
                MPb(li,lj) = -bb(li)*bb(lj) /mhdmuzero
              end if
            end do
          end do

          MPa(1,1) = dot_product(MPa(1,:), nwa(:)) + qa(1)
          MPa(2,1) = dot_product(MPa(2,:), nwa(:)) + qa(2)
          MPa(3,1) = dot_product(MPa(3,:), nwa(:)) + qa(3)

          MPb(1,1) = dot_product(MPb(1,:), nwb(:)) + qb(1)
          MPb(2,1) = dot_product(MPb(2,:), nwb(:)) + qb(2)
          MPb(3,1) = dot_product(MPb(3,:), nwb(:)) + qb(3)

          dva(:) = dva(:) - mb * ( &
                      MPa(:,1) / (rhoa**2 * oma) + &
                      MPb(:,1) / (rhob**2 * omb) &
                    ) - qd(:)

          dua   = dua + mb * ( &
                      dot_product(vab(:), nwa(:))*pa/(rhoa**2 * oma) - &
                      dot_product(vab(:), qa(:))/(rhoa * oma) + &
                      qc &
                    )

          dba(:) = dba(:) - 1./(rhoa * oma)*mb*( &
                      vab(:) * dot_product(ba(:),nwa(:)) - &
                      bb(:) * dot_product(vab(:),nwa(:)) &
                    )

        case (eeq_diffusion)
          if (difiso == 0) then
            do li = 1, 3
              do lj = 1,3
                ktb(li,lj) = difcond*bb(li)*bb(lj)
              end do
            end do
          end if

          call nw(rab, ra, rb, dr, ha, nwa)
          call nw(rab, ra, rb, dr, hb, nwb)

          if (s_ktp == esd_n2w) then
            call hessian(rab, ra, rb, ha, Hesa)
            do li = 1, 3
              do lj = 1,3
                tmpt1(li,lj) = tmpt1(li,lj) + mb/rhob * &
                                (ktb(li,lj) - kta(li,lj)) * nwa(li)
                tmpt2(li,lj) = tmpt2(li,lj) + mb/rhob * &
                                (tb - ta) * nwa(lj)
                tmpt3(li,lj) = tmpt3(li,lj) + kta(li,lj) *&
                                mb/rhob * (tb - ta) * Hesa(li,lj)
              end do
            end do
          else if ((s_ktp == esd_fw).or.(s_ktp == esd_fab)) then
            ktab(:,:) = (kta(:,:)+ktb(:,:))/2.
            call hessian(rab, ra, rb, ha, Hesa)
            ddta = ddta + mb/rhob * (tb - ta) * &
              ( dot_product(ktab(1,:),Hesa(1,:)) + &
                dot_product(ktab(2,:),Hesa(2,:)) + &
                dot_product(ktab(3,:),Hesa(3,:)) )
          else if (s_ktp == esd_2nw) then
            odda = 1./oma/rhoa/rhoa
            oddb = 1./omb/rhob/rhob
            qa(1) = dot_product(kta(1,:),dcftmp(:,i))
            qa(2) = dot_product(kta(2,:),dcftmp(:,i))
            qa(3) = dot_product(kta(3,:),dcftmp(:,i))
            qb(1) = dot_product(ktb(1,:),dcftmp(:,rj))
            qb(2) = dot_product(ktb(2,:),dcftmp(:,rj))
            qb(3) = dot_product(ktb(3,:),dcftmp(:,rj))

            ddta = ddta + mb*( &
              dot_product(qa(:),nwa(:))*odda + dot_product(qb(:),nwb(:))*oddb)
            else
              print*, "# <!> second deriv id not found"
              stop
          end if
        case (eeq_magnetohydrodiffusion)
          MPa(:,:) = 0.
          MPb(:,:) = 0.
          Mc(:)   = 0.
          odda = 1./oma/rhoa/rhoa
          oddb = 1./omb/rhob/rhob
          call nw(rab, ra, rb, dr, ha, nwa)
          call nw(rab, ra, rb, dr, hb, nwb)
          if (difiso == 0) then
            do li = 1, 3
              do lj = 1,3
                ktb(li,lj) = difcond*bb(li)*bb(lj)
              end do
            end do
          end if

          if (s_artts == 1) then
            call art_viscosity(rhoa, rhob, vab, urab, rab, dr, &
                                s_dim, ca, cb, ha, hb, qa, qb)
            call art_termcond(nwa, nwb, urab, pa, pb, ua, ub, rhoa, rhob, oma, omb, qc)
            call art_fdivbab(ba, bb, ma, nwa, nwb, oma, omb, rhoa, rhob, qd)
          end if

          Mc(1) = dot_product(ba(:),ba(:))
          Mc(2) = dot_product(bb(:),bb(:))
          do li = 1,3
            do lj = 1,3
              if (li == lj) then
                MPa(li,lj) = pa + (0.5*Mc(1) - ba(li)*ba(lj)) /mhdmuzero
                MPb(li,lj) = pb + (0.5*Mc(2) - bb(li)*bb(lj)) /mhdmuzero
              else
                MPa(li,lj) = -ba(li)*ba(lj) /mhdmuzero
                MPb(li,lj) = -bb(li)*bb(lj) /mhdmuzero
              end if
            end do
          end do

          MPa(1,1) = dot_product(MPa(1,:), nwa(:)) + qa(1)
          MPa(2,1) = dot_product(MPa(2,:), nwa(:)) + qa(2)
          MPa(3,1) = dot_product(MPa(3,:), nwa(:)) + qa(3)

          MPb(1,1) = dot_product(MPb(1,:), nwb(:)) + qb(1)
          MPb(2,1) = dot_product(MPb(2,:), nwb(:)) + qb(2)
          MPb(3,1) = dot_product(MPb(3,:), nwb(:)) + qb(3)

          dva(:) = dva(:) - mb * ( &
                      MPa(:,1) / (rhoa**2 * oma) + &
                      MPb(:,1) / (rhob**2 * omb) &
                    ) - qd(:)

          dua   = dua + mb * ( &
                      dot_product(vab(:), nwa(:))*pa/(rhoa**2 * oma) - &
                      dot_product(vab(:), qa(:))/(rhoa * oma) + &
                      qc &
                    )

          dba(:) = dba(:) - 1./(rhoa * oma)*mb*( &
                      vab(:) * dot_product(ba(:),nwa(:)) - &
                      bb(:) * dot_product(vab(:),nwa(:)) &
                    )
          ! diffusion
          if (s_ktp == esd_n2w) then
            call hessian(rab, ra, rb, ha, Hesa)
            do li = 1, 3
              do lj = 1,3
                tmpt1(li,lj) = tmpt1(li,lj) + mb/rhob * &
                                (ktb(li,lj) - kta(li,lj)) * nwa(li)
                tmpt2(li,lj) = tmpt2(li,lj) + mb/rhob * &
                                (tb - ta) * nwa(lj)
                tmpt3(li,lj) = tmpt3(li,lj) + kta(li,lj) *&
                                mb/rhob * (tb - ta) * Hesa(li,lj)
              end do
            end do
          else if ((s_ktp == esd_fw).or.(s_ktp == esd_fab)) then
            ktab(:,:) = (kta(:,:)+ktb(:,:))/2.
            call hessian(rab, ra, rb, ha, Hesa)
            ddta = ddta + mb/rhob * (tb - ta) * &
              ( dot_product(ktab(1,:),Hesa(1,:)) + &
                dot_product(ktab(2,:),Hesa(2,:)) + &
                dot_product(ktab(3,:),Hesa(3,:)) )
          else if (s_ktp == esd_2nw) then
            odda = 1./oma/rhoa/rhoa
            oddb = 1./omb/rhob/rhob
            qa(1) = dot_product(kta(1,:),dcftmp(:,i))
            qa(2) = dot_product(kta(2,:),dcftmp(:,i))
            qa(3) = dot_product(kta(3,:),dcftmp(:,i))
            qb(1) = dot_product(ktb(1,:),dcftmp(:,rj))
            qb(2) = dot_product(ktb(2,:),dcftmp(:,rj))
            qb(3) = dot_product(ktb(3,:),dcftmp(:,rj))

            ddta = ddta + mb*( &
              dot_product(qa(:),nwa(:))*odda + dot_product(qb(:),nwb(:))*oddb)
            else
              print*, "# <!> second deriv id not found"
              stop
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
        ! print*, 4, 95, dha + mb * dot_product(vab(:), nwa(:))
        if (s_adden == 1) then
          dha = dha + mb * dot_product(vab(:), nwa(:))
        end if
      end do overb

      if ( s_adden == 1 ) then
        dha =  (- ha / (s_dim * rhoa)) * dha / oma
      end if
      select case (s_ttp)
      case(eeq_hydro, eeq_magnetohydro)
      case(eeq_diffusion, eeq_magnetohydrodiffusion)
        if (s_ktp == esd_n2w) then
          do li = 1,3
            do lj = 1,3
              ddta = ddta + tmpt1(li,lj) * tmpt2(li,lj) + tmpt3(li,lj)
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
      store(es_ax:es_az,i) = dva(:)
      store(es_du,i) = dua
      store(es_dh,i) = dha
      store(es_ddt,i) = ddta
      store(es_dbx:es_dbz,i) = dba(:)
    end do overa
    !$omp end parallel do
    if (sum(store(es_ax:es_az,:)) /= 0.) then
      print*,sum(store(es_ax:es_az,:))
      read*
    end if
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
      urab, rab, dr, dim, ca, cb, ha, hb, &
      qa, qb)
    real, intent(in)  :: da, db, vab(3), urab(3), rab(3),&
                         dr, ca, cb, ha, hb
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
  !
  !   mrhoi5  = 0.5*pmassi*rho1i
  !
  !   runix = dx*rij1
  !   runiy = dy*rij1
  !   runiz = dz*rij1
  !   dvx = xpartveci(ivxi) - vxyzu(1,j)
  !   dvy = xpartveci(ivyi) - vxyzu(2,j)
  !   dvz = xpartveci(ivzi) - vxyzu(3,j)
  !   projv = dvx*runix + dvy*runiy + dvz*runiz
  !
  !   avBtermj = mrhoj5*alphaB*rho1j
  !   avBterm = mrhoi5*alphaB*rho1i
  !
  !   grkernj = grkern(q2j,qj)*hj21*hj21*cnormk*gradh(1,j)
  !   grkerni = grkern(q2i,qi)*hi41*cnormk*gradhi
  !   vsigB = sqrt((dvx - projv*runix)**2 + (dvy - projv*runiy)**2 + (dvz - projv*runiz)**2)
  !
  !   dBdissterm = (avBterm*grkerni + avBtermj*grkernj)*vsigB
  !
  !   dBx = Bxi - Bxj
  !   dBy = Byi - Byj
  !   dBz = Bzi - Bzj
  !   ! find Bevol and alphaB
  !   dB2 = dBx*dBx + dBy*dBy + dBz*dBz
  !   dudtresist = -0.5*dB2*dBdissterm
  ! end subroutine
end module
