module circuit2
  use omp_lib
  use timing,           only: addTime
  use kernel,           only: get_hessian, &
                              get_krad, &
                              get_n2w, &
                              get_nw
  use state,            only: get_difftype, &
                              getdim, &
                              get_tasktype, &
                              get_kerntype
  use BC
  use neighboursearch,  only: getneighbours,&
                              getNeibListL1,&
                              getNeibListL2

  implicit none

  public :: c2

  private
  save
    integer(8) :: start=0, finish=0

contains

  subroutine c2(c, ptype, pos, v, dv, mas, den, h, om, P, u, du, dh, cf, dcf, kcf, dfdx)
    real, allocatable, intent(in)    :: pos(:,:), v(:,:), mas(:), h(:), den(:), P(:), c(:),&
                                        u(:), cf(:), kcf(:), om(:)
    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(inout) :: dv(:,:), du(:), dh(:), dcf(:), dfdx(:,:,:)
    real                 :: dr, rhoa, rhob, qa, qb, qc, n2wa, n2wb, kr, r2, &
                            nwa(3), nwb(3), rab(3), vab(3), vba(3), urab(3), Pa(3), Pb(3), &
                            projv, df, ddf, Hes(3,3), oddi, oddj, tau2, r(3)
    integer, allocatable :: nlista(:), nlistb(:)
    integer              :: i, j, la, lb, n, dim, ttp, ktp, dtp, usekorrection
    integer(8)           :: t0, tneib


    ! usekorrection = 1
    usekorrection = 0
    call system_clock(start)
    n = size(ptype)
    tneib = 0.

    call getdim(dim)
    call get_krad(kr)
    call get_tasktype(ttp)
    call get_kerntype(ktp)
    call get_difftype(dtp)

    if (( ktp == 3 ).and.( dtp == 1 )) then
      call gradf(dim, ptype, mas, den, pos, v, h, om, dfdx)
    end if
    call getNeibListL1(nlista)

    !$omp parallel do default(none)&
    !$omp private(rab, dr, vab, urab, rhoa, rhob, nwa, nwb, qa, qb, qc, Pa, Pb)&
    !$omp private(n2wa, n2wb, j, i, r2, oddi ,oddj, la, lb)&
    !$omp private(projv, df, ddf, nlistb, Hes, vba, t0, tau2) &
    !$omp shared(dv, du, dh, dcf, n, pos, h, v, den, c, p, om, mas, u, kcf, cf)&
    !$omp shared(dim, kr, ktp, dtp, ttp, ptype, dfdx, nlista, usekorrection)&
    !$omp reduction(+:tneib)
    do la = 1, size(nlista)
      i = nlista(la)
      tau2 = 0.
      dv(:,i) = 0.
      du(i) = 0.
      dh(i) = 0.
      dcf(i) = 0.
      ! print*, i
      ! print*, pos(:,i)
      ! print*, dfdx(1,:,i)
      ! print*, dfdx(2,:,i)
      ! print*, dfdx(3,:,i)
      ! read*
      call getneighbours(i, pos, h, nlistb, t0)
      tneib = tneib + t0
      do lb = 1, size(nlistb)
        ! print *, i, nlist
        ! read *
        j = nlistb(lb)
        rab(:) = pos(:,i) - pos(:,j)
        r2 = dot_product(rab(:),rab(:))
        dr = sqrt(r2)
        vab(:) = v(:,i) - v(:,j)
        vba(:) = v(:,j) - v(:,i)
        urab(:) = rab(:) / dr
        select case (ttp)
        case (1, 9)
          ! hydroshock ! soundwave
          qa = 0.
          qb = 0.
          qc = 0.
          rhoa = den(i)
          rhob = den(j)
          ! print*, 100
          ! print*, i
          ! print*, rab
          ! print*, h(i)
          call get_nw(rab, h(i), nwa)
          ! print*, 999
          call get_nw(rab, h(j), nwb)
          ! call get_n2w(r, h(i), n2w)
          call art_viscosity(rhoa, rhob, vab, urab, c(i), c(j), qa, qb)
          call art_termcond(P(i), P(j), rhoa, rhob, qc)
          Pa(:) = (P(i) + qa) * nwa(:) / (rhoa**2 * om(i))
          Pb(:) = (P(j) + qb) * nwb(:) / (rhob**2 * om(j))

          dv(:,i) = dv(:,i) - mas(j) * (Pa(:) + Pb(:))

          du(i)   = du(i) + mas(j) * dot_product(vab(:),Pa(:)) &
                        + mas(j) / (0.5 *(rhoa + rhob)) * qc * (u(i) - u(j)) * &
                        0.5 * dot_product((nwa(:) + nwb(:)),urab(:))

          dh(i)   = dh(i) + mas(j) * dot_product(vab(:), nwa(:))
        case (2, 3)
          ! heatconduction
          call get_n2w(rab, h(i), n2wa)

          du(i) = du(i) - mas(j) / (den(i) * den(j)) * 2. * kcf(i) * kcf(j) &
                  / (kcf(i) + kcf(j)) * (cf(i) - cf(j)) * n2wa
          if (usekorrection == 1) then
            call LaplaceCorrection(rab, mas(j), den(j), dim, n2wa, tau2)
          end if
          ! du(i) = du(i) - mas(j) / (den(i) * den(j)) * (kcf(i) + kcf(j)) / 2. &
          !               * (cf(i) - cf(j)) * n2w
        case(4)
          ! photoevaporation
          ![ cf ~ eps ]![ dcf ~ deps/dt ]![ kcf ~ t_s]

          qa = 0.
          qb = 0.
          qc = 0.
          rhoa = den(i)
          rhob = den(j)

          call get_nw(rab, h(i), nwa)
          call get_nw(rab, h(j), nwb)
          call art_viscosity(rhoa, rhob, vab, urab, c(i), c(j), qa, qb)

          qa = qa * (1 - cf(i))
          qb = qb * (1 - cf(j))
          qc = 1/(1 - cf(i)) * mas(j) *&
          ( &
             + 0.25 * rhoa * abs(dot_product(vab(:), urab(:))) * (u(i) - u(j))/om(i)/rhoa/rhoa &
             * dot_product(urab,nwa)/dot_product(urab,urab) &
             + 0.25 * rhob * abs(dot_product(vab(:), urab(:))) * (u(j) - u(i))/om(j)/rhob/rhob &
             * dot_product(urab,nwb)/dot_product(urab,urab) &
          )
          ! dfgrhs is alerady sum, so it's not needed to sum it again.
          ! dv(:,i) = dfgrhs(:,i) * (1 + 1/(1 - cf(i)))
          !
          ! dcf(i)  = dcf(i) - mas(j) * &
          !         ( &
          !           cf(i)*(1 - cf(i))*kcf(i)/om(i)/rhoa*dot_product(-dfgrhs(:,i)/(1 - cf(i)),nwa(:)) + &
          !           cf(j)*(1 - cf(j))*kcf(j)/om(j)/rhob*dot_product(-dfgrhs(:,j)/(1 - cf(j)),nwb(:)) &
          !         )
          !
          ! du(i)   = du(i) + 1/om(i)/(1 - cf(i))/rhoa/rhoa * mas(j) * (P(i) + qa) * dot_product(vab(:),nwa(:)) -&
          !           cf(i) * kcf(i) / om(i) / rhoa * dot_product(-dfgrhs(:,i)/(1 - cf(i)), &
          !           mas(j) * (u(i) - u(j)) * nwa &
          !         + qc &
          !         )
          dh(i)   = dh(i) + mas(j) * dot_product(vab(:), nwa(:))
        case(5)
          ! 'diff-laplace'
          if (dtp == 1) then
            ! diff form
            if ( ktp /= 3 ) then
              ! n2w fab
              call get_n2w(rab, h(i), n2wa)
              if (usekorrection == 1) then
                call LaplaceCorrection(rab, mas(j), den(j), dim, n2wa, tau2)
              end if
              dv(:,i)  = dv(:,i) + mas(j)/den(j) * vba(:) * n2wa
            else
              ! 2nw
              call get_nw(rab, h(i), nwa)
              dv(1,i) = dv(1,i) + mas(j)/den(j) * ((dfdx(1,1,j) - dfdx(1,1,i))*nwa(1))
              dv(2,i) = dv(2,i) + mas(j)/den(j) * ((dfdx(2,2,j) - dfdx(2,2,i))*nwa(2))
              dv(3,i) = dv(3,i) + mas(j)/den(j) * ((dfdx(3,3,j) - dfdx(3,3,i))*nwa(3))
            end if
          elseif (dtp == 2) then
            ! symm form
            if ( ktp /= 3 ) then
              ! n2w fab
            else
              ! 2nw
              call get_nw(rab, h(i), nwa)
              call get_nw(rab, h(j), nwb)
              oddi = 1./om(i)/den(i)/den(i)
              oddj = 1./om(j)/den(j)/den(j)
              dv(1,i) = dv(1,i) + mas(j) * (dfdx(1,1,i)*nwa(1)*oddi + &
                                            dfdx(1,1,j)*nwb(1)*oddj)
              dv(2,i) = dv(2,i) + mas(j) * (dfdx(2,2,i)*nwa(2)*oddi + &
                                            dfdx(2,2,j)*nwb(2)*oddj)
              dv(3,i) = dv(3,i) + mas(j) * (dfdx(3,3,i)*nwa(3)*oddi + &
                                            dfdx(3,3,j)*nwb(3)*oddj)
            end if
          else
            print *, 'Diff type is not set in circuit2'
            stop
          end if
        case(6)
          ! diff-graddiv
          if (dtp == 1) then
            ! diff form
            if (ktp /= 3) then
              ! n2w fab
              call get_hessian(rab, h(i), Hes)
              dv(1,i) = dv(1,i) + mas(j)/den(j) * (vba(1)*Hes(1,1) + vba(2)*Hes(1,2) + vba(3)*Hes(1,3))
              dv(2,i) = dv(2,i) + mas(j)/den(j) * (vba(1)*Hes(2,1) + vba(2)*Hes(2,2) + vba(3)*Hes(2,3))
              dv(3,i) = dv(3,i) + mas(j)/den(j) * (vba(1)*Hes(3,1) + vba(2)*Hes(3,2) + vba(3)*Hes(3,3))
            else
              ! 2nw
              call get_nw(rab, h(i), nwa)
              ! qa = dfdx(1,1,i) + dfdx(2,2,i) + dfdx(3,3,i)
              ! qb = dfdx(1,1,j) + dfdx(2,2,j) + dfdx(3,3,j)
              qa = v(1,i)
              qb = v(1,j)
              dv(:,i) = dv(:,i) + mas(j) * (qb - qa) * nwa(:)
            end if
          elseif (dtp == 2) then
            ! symm form
            if ( ktp /= 3 ) then
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
        case default
          print *, 'Task type was not defined in circuit2.f90: line 240.'
          stop
        end select
      end do
      select case (ttp)
      case(1,2,3,4)
        dh(i) =  (- h(i) / (dim * den(i))) * dh(i) / om(i)
        if ( usekorrection == 1 ) then
          tau2 = tau2 / dim
          ! print*, tau2
          du(i) = du(i) / tau2
        end if
      case(5,6)
        ! diff-graddiv ! diff-laplace
        if ( usekorrection == 1 ) then
          tau2 = tau2 / dim
          print*, tau2
          dv(:,i) = dv(:,i) / tau2
        end if
        if ( ktp == 3 ) then
          dv(:,i) = dv(:,i) * den(i)
          ! print*, dv(:,i), den(i)
        end if
      case default
        print *, 'Task type was not set in circuit2 outside circle'
        stop
      end select
      ! print*, i, dv(:,i)
      ! read*
      ! print*, i, den(i)
    end do
    !$omp end parallel do
    ! print*, du(1:5)
    call system_clock(finish)
    call addTime(' circuit2', finish - start - tneib)
  end subroutine c2

  subroutine art_termcond(pa, pb, da, db, vsigu)
    real, intent(in)  :: pa, pb, da, db
    real, intent(out) :: vsigu

    vsigu = sqrt(abs(pa - pb)/(0.5 * (da + db)))
  end subroutine art_termcond

  subroutine art_viscosity(da, db, vab, urab, ca, cb, qa, qb)
    real, intent(in)  :: da, db, vab(3), urab(3), ca, cb
    real, intent(out) :: qa, qb
    real              :: alpha, betta, dvr
    qa = 0.
    qb = 0.
    alpha = 1.
    betta = 2.

    dvr = dot_product(vab,urab)
    if ( dvr < 0) then
      qa = -0.5 * da * (alpha*ca - betta*dvr) * dvr
      qb = -0.5 * db * (alpha*cb - betta*dvr) * dvr
    end if
  end subroutine art_viscosity

  subroutine gradf(dim, t, m, d, x, v, h, om, nv)
    real, allocatable, intent(in)  :: m(:), d(:), x(:,:), v(:,:), h(:), om(:)
    integer, allocatable, intent(in) :: t(:)
    real, allocatable, intent(inout) :: nv(:,:,:)
    integer, intent(in)            :: dim

    integer, allocatable :: nlista(:), nlistb(:)
    integer              :: n, i, j, la, lb, ni, nj, li
    integer(8)           :: t0, tneib

    real    :: vba(3), nw(3), rab(3), nwi(3), nwj(3), oddi, oddj

    n = size(m)
    ! if ( .not.allocated(nv) ) then
    !   allocate(nv(3,3,n))
    ! end if
    tneib = 0.
    ! call getNeibListL2(nlista)
    call getNeibListL1(nlista)

    !$omp parallel do default(none)&
    !$omp private(rab, vba, nw, i, j, la, lb, ni, nj, li, nlistb, t0) &
    !$omp private(nwi, nwj, oddi, oddj)&
    !$omp shared(dim, t, m, d, x, v, h, om, nv, n, nlista)&
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
  end subroutine gradf

  subroutine LaplaceCorrection(r, m, d, dim, n2w, tau2)
    real, intent(in)    :: r(3), m, d, n2w
    integer, intent(in) :: dim
    real, intent(inout) :: tau2

    tau2 = tau2 + 0.5 * m/d * r(1)*r(1) * n2w ! Hes(1,1)
    if ( dim > 1 ) then
      tau2 = tau2 + 0.5 * m/d * r(1)*r(2) * n2w ! Hes(1,2)
      tau2 = tau2 + 0.5 * m/d * r(2)*r(2) * n2w ! Hes(2,2)
      tau2 = tau2 + 0.5 * m/d * r(1)*r(2) * n2w ! Hes(1,2)
    end if
    if ( dim == 3 ) then
      tau2 = tau2 + 0.5 * m/d * r(1)*r(3) * n2w ! Hes(1,3)
      tau2 = tau2 + 0.5 * m/d * r(2)*r(3) * n2w ! Hes(2,3)
      tau2 = tau2 + 0.5 * m/d * r(3)*r(3) * n2w ! Hes(3,3)
      tau2 = tau2 + 0.5 * m/d * r(2)*r(3) * n2w ! Hes(2,3)
      tau2 = tau2 + 0.5 * m/d * r(1)*r(3) * n2w ! Hes(1,3)
    end if
  end subroutine
end module circuit2
