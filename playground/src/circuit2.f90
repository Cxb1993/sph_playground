module circuit2
  use const
  use omp_lib
  use timing,           only: addTime
  use kernel,           only: hessian, &
                              get_krad, &
                              n2w, &
                              nw, &
                              hessian_rr,&
                              get_w
  use neighboursearch,  only: getneighbours,&
                              getNeibListL1,&
                              getNeibListL2
  use state,            only: getdim, &
                              get_equations, getddwtype, &
                              getAdvancedDensity, &
                              getArtificialTerms, &
                              ginitvar, &
                              getdiffisotropic, &
                              getdiffconductivity, &
                              getmhdmagneticpressure,&
                              getPartNumber, &
                              getEqComponent
  use BC,               only: getCrossRef
  implicit none

  public :: c2init, c2

  private
  save
    integer(8)  :: start=0, finish=0
    integer :: &
      s_dim, s_ttp, s_ktp, s_adden, s_artts, s_ivt, initdone = 0,&
      eqSet(eqs_total)
    real        :: s_kr

contains

  subroutine c2init()
    call getdim(s_dim)
    call get_krad(s_kr)
    call get_equations(s_ttp)
    call getddwtype(s_ktp)
    call getAdvancedDensity(s_adden)
    call getArtificialTerms(s_artts)
    call ginitvar(s_ivt)
    call getEqComponent(eqSet)
    initdone = 1
  end subroutine

  subroutine c2(store, maxconsenrg)
    real, allocatable, intent(inout) :: store(:,:)
    real, intent(out) :: maxconsenrg

    real :: &
      qa(3), qb(3), qc, qd(3), qe, &
      dr, rhoa, rhob, ra(3), rb(3), r2, &
      rab(3), vab(3), vba(3), urab(3), &
      odda, oddb, kta(3,3), ktb(3,3), ktab(3,3), &
      tmpt1(3,3), tmpt2(3,3), tmpt3(3,3), MPa(3,3), MPb(3,3), Mc(3), &
      dva(3), dua, dha, ddta, dba(3), ba(3), bb(3), va(3), vb(3), &
      ha, hb, ca, cb, pa, pb, ua, ub, oma, omb, ta, tb, ma, mb, &
      wa, wb, nwa(3), nwb(3), n2wa, n2wb, Hesa(3,3), Hesb(3,3), &
      uba(3), ubb(3), dvaterm(3), dvbterm(3), &
      dtadx(3), dtbdx(3), kdtadx(3), kdtbdx(3), &
      difcond, mhdmuzero, &
      drhoadt, consenrg, &
      prada, pradb, radinteration

    integer, allocatable :: nlista(:), nlistb(:)
    integer :: &
      i, j, rj, la, lb, li, lj, difiso
    integer(8) :: t0, tneib

    call system_clock(start)
    call getmhdmagneticpressure(mhdmuzero)
    call getdiffisotropic(difiso)
    call getdiffconductivity(difcond)
    ! print*, eqSet
    ! read*

    if (initdone == 0) then
      call c2init()
    end if

    maxconsenrg = 0
    tneib = 0
    t0 = 0
    consenrg = 0.
    call getNeibListL1(nlista)
    !$omp parallel do default(none)&
    !$omp private(qa, qb, qc, qd, qe)&
    !$omp private(rab, dr, ra, rb, vab, urab, rhoa, rhob)&
    !$omp private(i, j, rj, r2, odda ,oddb, la, lb, drhoadt)&
    !$omp private(nlistb, vba, t0, kta, ktb, ktab, kdtadx, kdtbdx)&
    !$omp private(tmpt1, tmpt2, tmpt3, li, lj, MPa, MPb, Mc)&
    !$omp private(ha, hb, ca, cb, pa, pb, ua, ub, oma, omb, ta, tb, ma, mb)&
    !$omp private(wa, wb, n2wa, n2wb, nwa, nwb, Hesa, Hesb)&
    !$omp private(dva, dua, dha, ddta, dba, ba, bb, va, vb, dtadx, dtbdx)&
    !$omp private(uba, ubb, dvaterm, dvbterm)&
    !$omp private(prada, pradb, radinteration)&
    !$omp shared(store, nlista)&
    !$omp shared(eqSet)&
    !$omp shared(s_dim, s_kr, s_ktp, s_ttp, s_adden, s_artts, s_ivt)&
    !$omp shared(nw, hessian, hessian_rr)&
    !$omp shared(difcond, difiso, mhdmuzero)&
    !$omp reduction(+:tneib, consenrg)
    overa: do la = 1, size(nlista)
      i = nlista(la)
      rhoa  = store(es_den,i)
      ra(:) = store(es_rx:es_rz,i)
      ba(:) = store(es_bx:es_bz,i)
      uba(:) = ba(:)/sqrt(dot_product(ba(:),ba(:)))
      va(:) = store(es_vx:es_vz,i)
      ha  = store(es_h,i)
      ca  = store(es_c,i)
      pa  = store(es_p,i)
      ua  = store(es_u,i)
      oma = store(es_om,i)
      ta  = store(es_t,i)
      ma  = store(es_m,i)
      dtadx(:) = store(es_dtdx:es_dtdz,i)
      kdtadx(:) = 0.
      dva(:) = 0.
      dua = 0.
      dha = 0.
      ddta = 0.
      dba(:) = 0.
      drhoadt = 0.

      tmpt1(:,:) = 0.
      tmpt2(:,:) = 0.
      tmpt3(:,:) = 0.
      kta(:,:) = 0.
      ktb(:,:) = 0.

      MPa(:,:) = 0.
      Mc(1) = dot_product(ba(:),ba(:))

      prada = 0.
      pradb = 0.
      radinteration = 0.

      if (eqSet(eqs_magneto) == 1) then
        do li = 1,3
          MPa(li,li) = (0.5*Mc(1) - ba(li)*ba(li))/mhdmuzero
          do lj = 1,3
            if (li /= lj) then
              MPa(li,lj) = -ba(li)*ba(lj)/mhdmuzero
            end if
          end do
        end do
        ! diffusion done with respect to magnetic fields
        if (eqSet(eqs_diff) == 1) then
          if (difiso == 1) then
            kta(1,1) = difcond
            kta(2,2) = difcond
            kta(3,3) = difcond

            ktb(1,1) = difcond
            ktb(2,2) = difcond
            ktb(3,3) = difcond
          else
            ! mti-like diffusion only for middle layer
            ! if ((ra(2) >= 0.).and.(ra(2) < 1./3.)) then
            !   kta(1,1) = difcond
            !   kta(2,2) = difcond
            !   kta(3,3) = difcond
            ! else if (ra(2) > 1.) then
            !   kta(1,1) = difcond
            !   kta(2,2) = difcond
            !   kta(3,3) = difcond
            ! else
            !   do li = 1,3
            !     do lj = 1,3
            !       kta(li,lj) = difcond*uba(li)*uba(lj)
            !     end do
            !   end do
            ! end if
            do li = 1, 3
                kta(li,li) = difcond*uba(li)*uba(li)
            end do
          end if
        end if
      end if

      if (eqSet(eqs_diff) == 1) then
        ! if (doLimiter == 1) then
        !   ! R = abs(grad(rho*ksi))/kappa/rho^2/ksi
        !   fluxlimR = 1.
        !   ! state(es_fluxlim,i) = (2. + R)/(6. + 3.*R + R*R)
        !   state(es_fluxlim,i) = 1.
        ! end if

        if (difiso == 1) then
          kta(1,1) = difcond
          kta(2,2) = difcond
          kta(3,3) = difcond

          ktb(1,1) = difcond
          ktb(2,2) = difcond
          ktb(3,3) = difcond
        ! else
        !   do li = 1, 3
        !       kta(li,li) = difcond*uba(li)*uba(li)
        !   end do
        end if
      end if

      if (eqSet(eqs_fld) == 1) then
        prada = 1./3.*rhoa*ta
        radinteration = c * kappa * rhoa * ta - 4. * kappa * sigma_b * (u(i) / cv)**4
        dua  = radinteration
        ddta = -radinteration*rhoa
      end if

      call getneighbours(i, nlistb, t0)
      tneib = tneib + t0
      overb: do lb = 1, size(nlistb)
        j = nlistb(lb)
        rj = getCrossRef(j)

        rb(:) = store(es_rx:es_rz, j)
        bb(:) = store(es_bx:es_bz, rj)
        ubb(:) = bb(:)/sqrt(dot_product(bb(:),bb(:)))
        vb(:) = store(es_vx:es_vz, rj)
        rhob = store(es_den, rj)
        hb = store(es_h, rj)
        cb = store(es_c, rj)
        pb = store(es_p, rj)
        ub = store(es_u, rj)
        omb = store(es_om, rj)
        tb = store(es_t, rj)
        mb = store(es_m, rj)
        dtbdx(:) = store(es_dtdx:es_dtdz, rj)
        kdtbdx(:) = 0.
        qa(:) = 0.
        qb(:) = 0.
        qc    = 0.
        qd(:) = 0.
        qe    = 0.
        rab(:) = ra(:) - rb(:)
        r2 = dot_product(rab(:),rab(:))
        dr = sqrt(r2)
        vab(:) = va(:) - vb(:)
        vba(:) = vb(:) - va(:)
        urab(:) = rab(:)/dr
        Hesb(:,:) = 0.

        call nw(rab, ra(:), rb(:), dr, ha, nwa)
        call nw(rab, ra(:), rb(:), dr, hb, nwb)
        ! print*, nwa
        ! print*, nwb
        ! call get_w(dr, ha, wa)
        ! call get_w(dr, hb, wb)
        ! nwa(:) = -2.*(239./231.)*wa*rab(:)/ha/ha
        ! nwb(:) = -2.*(239./231.)*wb*rab(:)/hb/hb
        ! nwa(:) = -2.*(1.5)*wa*rab(:)/ha/ha
        ! nwb(:) = -2.*(1.5)*wb*rab(:)/hb/hb
        ! print*, nwa
        ! print*, nwb
        ! read*
        drhoadt = drhoadt + mb * dot_product(vab(:),nwa(:))
        if (s_adden == 1) then
          dha = dha + mb * dot_product(vab(:), nwa(:))
        end if

        if (eqSet(eqs_fld) == 1) then
          pradb = 1./3.*rhob*tb

          if (s_ktp /= esd_2nw) then
            ktab(:,:) = (kta(:,:)+ktb(:,:))/2.
            ! call hessian(rab, ra, rb, ha, Hesa)
            call n2w(rab, ha, n2wa)
            call n2w(rab, hb, n2wb)
            ddta = ddta + mb/rhob/2. * (Da * n2wa + Db * n2wb) * (prada - pradb) &
                        - prada/rhoa*dot_product(vba,nwa)
          else
            odda = 1./oma/rhoa/rhoa
            oddb = 1./omb/rhob/rhob
            kdtadx(1) = dot_product(kta(1,:),dtadx(:))
            kdtadx(2) = dot_product(kta(2,:),dtadx(:))
            kdtadx(3) = dot_product(kta(3,:),dtadx(:))
            kdtbdx(1) = dot_product(ktb(1,:),dtbdx(:))
            kdtbdx(2) = dot_product(ktb(2,:),dtbdx(:))
            kdtbdx(3) = dot_product(ktb(3,:),dtbdx(:))

            ddta = ddta + mb*( &
              dot_product(kdtadx(:),nwa(:))*odda + dot_product(kdtbdx(:),nwb(:))*oddb)
          end if
        end if

        if (eqSet(eqs_hydro) == 1) then
          if (s_artts == 1) then
            call art_viscosity(rhoa, rhob, vab, urab, rab, dr, s_dim, ca, cb, ha, hb, qa, qb)
            call art_termcond(nwa, nwb, urab, pa, pb, ua, ub, rhoa, rhob, oma, omb, qc)
          end if

          dva(:) = dva(:) - mb * (&
                      ((pa + prada) * nwa(:) + qa(:)) / (rhoa**2 * oma) + &
                      ((pb + pradb) * nwb(:) + qb(:)) / (rhob**2 * omb) &
                    )

          dua   = dua + mb * (&
                      ( dot_product(vab(:), pa * nwa(:)) +&
                        dot_product(vab(:), qa(:)) )&
                      /(rhoa**2 * oma) +&
                      qc &
                    )
        end if

        if (eqSet(eqs_diff) == 1) then
          if (difiso == 0) then
            ! mti-like diffusion definiotion
            ! if ((ra(2) > 0.).and.(ra(2) < 1./3.)) then
            !   ktb(1,1) = difcond
            !   ktb(2,2) = difcond
            !   ktb(3,3) = difcond
            ! else if (ra(2) > 1.) then
            !   ktb(1,1) = difcond
            !   ktb(2,2) = difcond
            !   ktb(3,3) = difcond
            ! else
            !   do li = 1,3
            !     do lj = 1,3
            !       ktb(li,lj) = difcond*ubb(li)*ubb(lj)
            !     end do
            !   end do
            ! end if
            ! regular anisodiffusion
            do li = 1, 3
                ktb(li,li) = difcond*ubb(li)*ubb(li)
            end do
          end if
          if (s_ktp /= esd_2nw) then
          ! if (s_ktp == esd_n2w) then
          !   call hessian(rab, ra, rb, ha, Hesa)
          !   do li = 1, 3
          !     do lj = 1, 3
          !       tmpt1(li,lj) = tmpt1(li,lj) + mb/rhob * &
          !                       (ktb(li,lj) - kta(li,lj)) * nwa(li)
          !       tmpt2(li,lj) = tmpt2(li,lj) + mb/rhob * &
          !                       (tb - ta) * nwa(lj)
          !       tmpt3(li,lj) = tmpt3(li,lj) + kta(li,lj) * &
          !                       mb/rhob * (tb - ta) * Hesa(li,lj)
          !     end do
          !   end do
          ! else if ((s_ktp == esd_fw).or.(s_ktp == esd_fab)) then
            ktab(:,:) = (kta(:,:)+ktb(:,:))/2.
            call hessian(rab, ra, rb, ha, Hesa)
            ! call hessian(rab, ra, rb, hb, Hesb)
            ! Hesa(:,:) = 0.5*(Hesa(:,:)+Hesb(:,:))
            ddta = ddta + mb/rhob * (tb - ta) * &
               (dot_product(ktab(1,:),Hesa(1,:)) + &
                dot_product(ktab(2,:),Hesa(2,:)) + &
                dot_product(ktab(3,:),Hesa(3,:)))
          else if (s_ktp == esd_2nw) then
            odda = 1./oma/rhoa/rhoa
            oddb = 1./omb/rhob/rhob
            kdtadx(1) = dot_product(kta(1,:),dtadx(:))
            kdtadx(2) = dot_product(kta(2,:),dtadx(:))
            kdtadx(3) = dot_product(kta(3,:),dtadx(:))
            kdtbdx(1) = dot_product(ktb(1,:),dtbdx(:))
            kdtbdx(2) = dot_product(ktb(2,:),dtbdx(:))
            kdtbdx(3) = dot_product(ktb(3,:),dtbdx(:))

            ddta = ddta + mb*( &
              dot_product(kdtadx(:),nwa(:))*odda + dot_product(kdtbdx(:),nwb(:))*oddb)
            else
              print*, "# <!> second deriv id not found"
              stop
          end if
        end if

        if (eqSet(eqs_magneto) == 1) then
          MPb(:,:) = 0.
          Mc(:)    = 0.
          odda = 1./oma/rhoa/rhoa
          oddb = 1./omb/rhob/rhob

          if (s_artts == 1) then
            call art_fdivbab(ba, bb, ma, nwa, nwb, oma, omb, rhoa, rhob, qd)
            call art_resistivity(nwa, nwb, Ba, Bb, vab, urab, rhoa, rhob, oma, omb, qe)
          end if

          Mc(2) = dot_product(bb(:),bb(:))
          do li = 1,3
            MPb(li,li) = (0.5*Mc(2) - bb(li)*bb(li))/mhdmuzero
            do lj = 1,3
              if (li /= lj) then
                MPb(li,lj) = -bb(li)*bb(lj)/mhdmuzero
              end if
            end do
          end do

          dvaterm(1) = dot_product(MPa(1,:), nwa(:))
          dvaterm(2) = dot_product(MPa(2,:), nwa(:))
          dvaterm(3) = dot_product(MPa(3,:), nwa(:))

          dvbterm(1) = dot_product(MPb(1,:), nwb(:))
          dvbterm(2) = dot_product(MPb(2,:), nwb(:))
          dvbterm(3) = dot_product(MPb(3,:), nwb(:))

          dva(:) = dva(:) - mb * ( &
                      dvaterm(:) / (rhoa**2 * oma) + &
                      dvbterm(:) / (rhob**2 * omb) &
                    ) - qd(:)

          dua    = dua + mb * qe

          dba(:) = dba(:) - 1./(rhoa * oma)*mb*( &
                      vab(:) * dot_product(ba(:),nwa(:)) - &
                      bb(:) * dot_product(vab(:),nwa(:)) &
                    )
        end if
      end do overb

      if ( s_adden == 1 ) then
        dha =  (- ha / (s_dim * rhoa)) * dha / oma
      end if

      if (eqSet(eqs_diff) == 1) then
        ! if (s_ktp == esd_n2w) then
          do li = 1,3
            do lj = 1,3
              ddta = ddta + tmpt1(li,lj) * tmpt2(li,lj) + tmpt3(li,lj)
            end do
          end do
        ! end if
      end if

      store(es_ax:es_az,i) = dva(:)
      store(es_dh,i) = dha
      store(es_ddt,i) = ddta/rhoa
      store(es_dbx:es_dbz,i) = dba(:)
      store(es_du,i) = dua
      ! store(es_du,i) = dua + ddta/rhoa
      consenrg = consenrg + ma*(dot_product(va(:),dva(:)) + dua + &
        dot_product(ba(:),dba(:))/rhoa - &
        0.5*dot_product(ba(:),ba(:))/rhoa/rhoa*drhoadt/oma)
    end do overa
    !$omp end parallel do
    if (abs(consenrg) > abs(maxconsenrg)) maxconsenrg = consenrg
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

  pure subroutine art_resistivity(nwa, nwb, Ba, Bb, vab, urab, rhoa, rhob, oma, omb, qe)
    real, intent(in)  :: &
      nwa(3), nwb(3), Ba(3), Bb(3), vab(3), urab(3), &
      rhoa, rhob, oma, omb
    real, intent(out) :: qe
    real              :: vsigB, Bab(3), vabxurab(3)
    qe = 0.
    Bab(:) = Ba(:) - Bb(:)
    vabxurab(1) = vab(2)*urab(3) - vab(3)*urab(2)
    vabxurab(2) = vab(3)*urab(1) - vab(1)*urab(3)
    vabxurab(3) = vab(1)*urab(2) - vab(2)*urab(1)
    vsigB = sqrt(dot_product(vabxurab(:),vabxurab(:)))
    qe = -0.25*(vsigB*dot_product(urab(:),nwa(:))/rhoa/rhoa/oma + &
                vsigB*dot_product(urab(:),nwb(:))/rhob/rhob/omb&
              )*dot_product(Bab(:),Bab(:))
  end subroutine
end module
