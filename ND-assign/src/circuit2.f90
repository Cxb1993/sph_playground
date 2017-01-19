module circuit2
  use omp_lib
  use kernel
  use BC

  implicit none

  public :: c2

  private

contains

  subroutine c2(n, c, pos, v, dv, mas, den, h, om, P, u, du, dh, cf, dcf, kcf)
    integer, intent(in) :: n
    real, allocatable, intent(in) :: pos(:,:), v(:,:), mas(:), h(:), den(:), P(:), c(:), om(:),&
                                     u(:), cf(:), kcf(:)
    real, allocatable, intent(inout) :: dv(:,:), du(:), dh(:), dcf(:)
    real                :: dr, rhoa, rhob, qa, qb, qc, n2w, kr, r2
    real                :: nwa(3), nwb(3), r(3), vab(3), urab(3), Pa(3), Pb(3), dfgrhs(3,n)
    integer             :: i, j, dim
    integer             :: tt!, kt

    call get_dim(dim)
    call get_tasktype(tt)
    call get_krad(kr)

    if (tt == 4) then
      call diff_force(mas, den, pos, v, h, c, P, om, dfgrhs)
    end if

    !$omp parallel do default(none)&
    !$omp private(r, dr, vab, urab, rhoa, rhob, nwa, nwb, qa, qb, qc, Pa, Pb, n2w, j, i, r2) &
    !$omp shared(dv, du, dh, dcf, n, pos, h, tt, v, den, c, p, om, mas, u, kcf, cf)&
    !$omp shared(dfgrhs, dim, kr)
    do i = 1, n
      dv(:,i) = 0.
      du(i) = 0.
      dh(i) = 0.
      dcf(i) = 0.
      do j = 1, n
        if (i /= j) then
          r(:) = pos(:,i) - pos(:,j)
          r2 = dot_product(r(:),r(:))
          if ((r2 < (kr * h(i))**2).or.(r2 < (kr * h(j))**2)) then
            dr = sqrt(r2)

            select case (tt)
            case (1)
              qa = 0.
              qb = 0.
              qc = 0.
              vab(:) = v(:,i) - v(:,j)
              urab(:) = r(:) / dr
              rhoa = den(i)
              rhob = den(j)

              call get_nw(r, h(i), nwa)
              call get_nw(r, h(j), nwb)
              call get_n2w(r, h(i), n2w)

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
              call get_n2w(r, h(i), n2w)

              du(i) = du(i) - mas(j) / (den(i) * den(j)) * 2. * kcf(i) * kcf(j) &
                      / (kcf(i) + kcf(j)) * (cf(i) - cf(j)) * n2w
              ! du(i) = du(i) - mas(j) / (den(i) * den(j)) * (kcf(i) + kcf(j)) / 2. &
              !               * (cf(i) - cf(j)) * n2w
            case(4)
              ! photoevaporation
              ![ cf ~ eps ]![ dcf ~ deps/dt ]![ kcf ~ t_s]

              qa = 0.
              qb = 0.
              qc = 0.
              vab(:) = v(:,i) - v(:,j)
              urab(:) = r(:) / dr
              rhoa = den(i)
              rhob = den(j)

              call get_nw(r, h(i), nwa)
              call get_nw(r, h(j), nwb)
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
              dv(:,i) = dfgrhs(:,i) * (1 + 1/(1 - cf(i)))

              dcf(i)  = dcf(i) - mas(j) * &
                      ( &
                        cf(i)*(1 - cf(i))*kcf(i)/om(i)/rhoa*dot_product(-dfgrhs(:,i)/(1 - cf(i)),nwa(:)) + &
                        cf(j)*(1 - cf(j))*kcf(j)/om(j)/rhob*dot_product(-dfgrhs(:,j)/(1 - cf(j)),nwb(:)) &
                      )

              du(i)   = du(i) + 1/om(i)/(1 - cf(i))/rhoa/rhoa * mas(j) * (P(i) + qa) * dot_product(vab(:),nwa(:)) -&
                        cf(i) * kcf(i) / om(i) / rhoa * dot_product(-dfgrhs(:,i)/(1 - cf(i)), &
                        mas(j) * (u(i) - u(j)) * nwa &
                      + qc &
                      )

              dh(i)   = dh(i) + mas(j) * dot_product(vab(:), nwa(:))
            case(5)
              ! 'diff-laplace'
              call get_n2w(r, h(i), n2w)
              dv(1,i)  = dv(1,i) + mas(j)/den(j) * (v(1,j) - v(1,i)) * n2w
            case default
              print *, 'Task type was not defined in circuit2'
              stop
            end select
          end if
        end if
      end do
      if (tt /= 5) then
        dh(i) =  (- h(i) / (dim * den(i))) * dh(i) / om(i)
      end if
    end do
    !$omp end parallel do
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

  subroutine diff_force(mas, den, pos, v, h, c, P, om, dfgrhs)
    real, intent(in)    :: mas(:), den(:), pos(:,:), v(:,:), h(:), c(:), P(:), om(:)
    real, intent(inout) :: dfgrhs(:,:)
    real                :: r(3), vab(3), urab(3), nwa(3), nwb(3), Pa(3), Pb(3),&
                           r2, dr, qa, qb, qc, rhoa, rhob, kr
    integer             :: i, j, szd1, szd2

    szd1 = 3
    szd2 = size(dfgrhs, dim=2)

    call get_krad(kr)

    !$omp parallel do default(none)&
    !$omp private(i, j, r, r2, dr, qa, qb, qc, vab, urab, rhoa, rhob, nwa, nwb, Pa, Pb)&
    !$omp shared(kr, mas, den, pos, v, h, c, P, om, dfgrhs, szd2)
    do i=1,szd2
      dfgrhs(:,i) = 0.
      do j=1,szd2
        if (i /= j) then
          r(:) = pos(:,i) - pos(:,j)
          r2 = dot_product(r(:),r(:))
          if ((r2 < (kr * h(i))**2).or.(r2 < (kr * h(j))**2)) then
            dr = sqrt(r2)
            qa = 0.
            qb = 0.
            qc = 0.
            vab(:) = v(:,i) - v(:,j)
            urab(:) = r(:) / dr
            rhoa = den(i)
            rhob = den(j)

            call get_nw(r, h(i), nwa)
            call get_nw(r, h(j), nwb)

            call art_viscosity(rhoa, rhob, vab, urab, c(i), c(j), qa, qb)
            call art_termcond(P(i), P(j), rhoa, rhob, qc)
            Pa(:) = (P(i) + qa) * nwa(:) / (rhoa**2 * om(i))
            Pb(:) = (P(j) + qb) * nwb(:) / (rhob**2 * om(j))

            dfgrhs(:,i) = dfgrhs(:,i) - mas(j) * (Pa(:) + Pb(:))
          end if
        end if
      end do
    end do
    !$omp end parallel do
  end subroutine diff_force
end module circuit2
