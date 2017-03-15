module circuit2
  use omp_lib
  use kernel
  use BC
  use neighboursearch, only: getneighbours,&
                             getNeibListL2

  implicit none

  public :: c2, setStepsize, getTime

  private
    integer, save :: stepsize = 1
    real :: start=0., finish=0., elapsed=0.

contains
  subroutine setStepsize(i)
    integer, intent(in) :: i
    stepsize = i
  end subroutine setStepsize

  subroutine getTime(ot)
    real, intent(out) :: ot
    ot = elapsed
  end subroutine getTime

  subroutine c2(c, ptype, pos, v, dv, mas, den, h, om, P, u, du, dh, cf, dcf, kcf, dfdx)
    real, allocatable, intent(in)    :: pos(:,:), v(:,:), mas(:), h(:), den(:), P(:), c(:),&
                                        u(:), cf(:), kcf(:), om(:)
    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(inout) :: dv(:,:), du(:), dh(:), dcf(:), dfdx(:,:,:)
    real                 :: dr, rhoa, rhob, qa, qb, qc, n2wa, n2wb, kr, r2, &
                            nwa(3), nwb(3), rab(3), vab(3), vba(3), urab(3), Pa(3), Pb(3), &
                            projv, df, ddf, Hes(3,3), t0, tneib
    integer, allocatable :: nlist(:)
    integer              :: i, j, l, n, dim, ttp, ktp, dtp

    call cpu_time(start)
    n = size(ptype)
    tneib = 0.

    call get_dim(dim)
    call get_krad(kr)
    call get_tasktype(ttp)
    call get_kerntype(ktp)
    call get_difftype(dtp)

    if (( ktp == 3 ).and.( dtp == 1 )) then
      call gradf(dim, ptype, mas, den, pos, v, h, dfdx)
    end if

    !$omp parallel do default(none)&
    !$omp private(rab, dr, vab, urab, rhoa, rhob, nwa, nwb, qa, qb, qc, Pa, Pb)&
    !$omp private(n2wa, n2wb, j, i, r2)&
    !$omp private(projv, df, ddf, nlist, Hes, vba, t0) &
    !$omp shared(dv, du, dh, dcf, n, pos, h, v, den, c, p, om, mas, u, kcf, cf)&
    !$omp shared(dim, kr, ktp, dtp, ttp, ptype, stepsize, dfdx)&
    !$omp reduction(+:tneib)
    do i = 1, n, stepsize
      ! print *, i, stepsize, ptype(i)
      ! read *
      if ( ptype(i) /= 0 ) then
        dv(:,i) = 0.
        du(i) = 0.
        dh(i) = 0.
        dcf(i) = 0.
        call getneighbours(i, pos, h, nlist, t0)
        tneib = tneib + t0
        do l = 1, size(nlist)
          ! print *, i, nlist
          ! read *
          j = nlist(l)
          rab(:) = pos(:,i) - pos(:,j)
          r2 = dot_product(rab(:),rab(:))
          dr = sqrt(r2)
          vab(:) = v(:,i) - v(:,j)
          vba(:) = v(:,j) - v(:,i)
          urab(:) = rab(:) / dr
          select case (ttp)
          case (1)
            qa = 0.
            qb = 0.
            qc = 0.
            rhoa = den(i)
            rhob = den(j)

            call get_nw(rab, h(i), nwa)
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
            call get_n2w(rab, h(i), n2wa)

            du(i) = du(i) - mas(j) / (den(i) * den(j)) * 2. * kcf(i) * kcf(j) &
                    / (kcf(i) + kcf(j)) * (cf(i) - cf(j)) * n2wa
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
              if ( ktp /= 3 ) then
                call get_n2w(rab, h(i), n2wa)
                dv(:,i)  = dv(:,i) + mas(j)/den(j) * vba(:) * n2wa
              else
                call get_nw(rab, h(i), nwa)
                dv(1,i) = dv(1,i) + mas(j)/den(j) * 0.5*((2*dfdx(1,1,j)-2*dfdx(1,1,i))*nwa(1) + &
                              (dfdx(1,2,j) + dfdx(2,1,j) - dfdx(1,2,i) - dfdx(2,1,i))*nwa(2) + &
                              (dfdx(1,3,j) + dfdx(3,1,j) - dfdx(1,3,i) - dfdx(3,1,i))*nwa(3))
                dv(2,i) = dv(2,i) + mas(j)/den(j) * 0.5*((dfdx(2,1,j) + dfdx(1,2,j) - dfdx(2,1,i) - dfdx(1,2,i))*nwa(1) + &
                              (2*dfdx(2,2,j) - 2*dfdx(2,2,i))*nwa(2) + &
                              (dfdx(2,3,j) + dfdx(3,2,j) - dfdx(2,3,i) - dfdx(3,2,i))*nwa(3))
                dv(3,i) = dv(3,i) + mas(j)/den(j) * 0.5*((dfdx(3,1,j) + dfdx(1,3,j) - dfdx(3,1,i) - dfdx(1,3,i))*nwa(1) + &
                              (dfdx(3,2,j) + dfdx(2,3,j) - dfdx(3,2,i) - dfdx(2,3,i))*nwa(2) + &
                              (2*dfdx(3,3,j) - 2*dfdx(3,3,i))*nwa(3))
              end if
            elseif (dtp == 2) then
              print*, 'symm'
            else
              print *, 'Diff type is not set in circuit2'
              stop
            end if
          case(6)
            ! diff-graddiv
            if (dtp == 1) then
              if (ktp /= 3) then
                call get_hessian(rab, h(i), Hes)
                dv(1,i) = dv(1,i) + mas(j)/den(j) * (vba(1)*Hes(1,1) + vba(2)*Hes(1,2) + vba(3)*Hes(1,3))
                dv(2,i) = dv(2,i) + mas(j)/den(j) * (vba(1)*Hes(2,1) + vba(2)*Hes(2,2) + vba(3)*Hes(2,3))
                dv(3,i) = dv(3,i) + mas(j)/den(j) * (vba(1)*Hes(3,1) + vba(2)*Hes(3,2) + vba(3)*Hes(3,3))
              else
                call get_nw(rab, h(i), nwa)
                dv(1,i) = dv(1,i) + mas(j)/den(j) * 0.5*((2*dfdx(1,1,j)-2*dfdx(1,1,i))*nwa(1) + &
                              (dfdx(1,2,j) + dfdx(2,1,j) - dfdx(1,2,i) - dfdx(2,1,i))*nwa(2) + &
                              (dfdx(1,3,j) + dfdx(3,1,j) - dfdx(1,3,i) - dfdx(3,1,i))*nwa(3))
                dv(2,i) = dv(2,i) + mas(j)/den(j) * 0.5*((dfdx(2,1,j) + dfdx(1,2,j) - dfdx(2,1,i) - dfdx(1,2,i))*nwa(1) + &
                              (2*dfdx(2,2,j) - 2*dfdx(2,2,i))*nwa(2) + &
                              (dfdx(2,3,j) + dfdx(3,2,j) - dfdx(2,3,i) - dfdx(3,2,i))*nwa(3))
                dv(3,i) = dv(3,i) + mas(j)/den(j) * 0.5*((dfdx(3,1,j) + dfdx(1,3,j) - dfdx(3,1,i) - dfdx(1,3,i))*nwa(1) + &
                              (dfdx(3,2,j) + dfdx(2,3,j) - dfdx(3,2,i) - dfdx(2,3,i))*nwa(2) + &
                              (2*dfdx(3,3,j) - 2*dfdx(3,3,i))*nwa(3))
              end if
            elseif (dtp == 2) then
              print*, 'symm'
            else
              print *, 'Diff type is not set in circuit2'
              stop
            end if
          case default
            print *, 'Task type was not defined in circuit2 inside circle'
            stop
          end select
        end do
        select case (ttp)
        case(1,2,3,4)
          dh(i) =  (- h(i) / (dim * den(i))) * dh(i) / om(i)
        case(5,6)
        case default
          print *, 'Task type was not set in circuit2 outside circle'
          stop
        end select
      end if
      ! print*, i, dv(:,i)
      ! read*
    end do
    !$omp end parallel do
    call cpu_time(finish)
    elapsed = elapsed + (finish - start) - tneib
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

  subroutine gradf(dim, t, m, d, x, v, h, nv)
    real, allocatable, intent(in)  :: m(:), d(:), x(:,:), v(:,:), h(:)
    integer, allocatable, intent(in) :: t(:)
    real, allocatable, intent(out) :: nv(:,:,:)
    integer, intent(in)            :: dim

    integer, allocatable :: nlista(:), nlistb(:)
    integer :: n, i, j, la, lb, ni, nj, li
    real    :: vba(3), nw(3), rab(3), t0, tneib
    call cpu_time(start)

    n = size(m)
    allocate(nv(3,3,n))
    nv(:,:,:) = 0.
    tneib = 0.
    call getNeibListL2(nlista)

    !$omp parallel do default(none)&
    !$omp private(rab, vba, nw, i, j, la, lb, ni, nj, li, nlistb, t0) &
    !$omp shared(dim, t, m, d, x, v, h, nv, n, stepsize, nlista)&
    !$omp reduction(+:tneib)
    do la = 1,size(nlista)
      i = nlista(la)
      call getneighbours(i, x, h, nlistb, t0)
      tneib = tneib + t0
      do lb = 1, size(nlistb)
        j = nlistb(lb)
        vba(:) = v(:,j) - v(:,i)
        rab(:) = x(:,i) - x(:,j)
        call get_nw(rab, h(i), nw)
        do ni = 1,dim
          do nj = 1,dim
            nv(ni,nj,i) = nv(ni,nj,i) + m(j)/d(j)*vba(ni)*nw(nj)
          end do
        end do
      end do
    end do
    !$omp end parallel do
    call cpu_time(finish)
    elapsed = elapsed + (finish - start) - tneib
  end subroutine gradf
end module circuit2
