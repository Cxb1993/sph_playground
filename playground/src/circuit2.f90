module circuit2
  use omp_lib
  use kernel
  use BC
  use neighboursearch, only:getneighbours

  implicit none

  public :: c2, setStepsize

  private
  integer, save :: stepsize = 1

contains
  subroutine setStepsize(i)
    integer, intent(in) :: i
    stepsize = i
  end subroutine setStepsize

  subroutine c2(c, ptype, pos, v, dv, mas, den, h, om, P, u, du, dh, cf, dcf, kcf)
    real, allocatable, intent(in)    :: pos(:,:), v(:,:), mas(:), h(:), den(:), P(:), c(:), om(:),&
                                        u(:), cf(:), kcf(:)
    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(inout) :: dv(:,:), du(:), dh(:), dcf(:)
    real                 :: dr, rhoa, rhob, qa, qb, qc, n2wa, n2wb, kr, r2, &
                            nwa(3), nwb(3), rab(3), vab(3), vba(3), urab(3), Pa(3), Pb(3), &
                            projv, df, ddf, Jac(3,3), sij
    integer, allocatable :: nlist(:)
    real, allocatable    :: gradv(:,:,:)
    integer              :: i, j, l, n, dim
    integer              :: tt, kt

    n = size(ptype)

    call get_dim(dim)
    call get_tasktype(tt)
    call get_krad(kr)
    call get_kerntype(kt)

    if (kt == 3) then
      call fungradient(dim, ptype, mas, den, pos, v, h, gradv)
    end if

    !$omp parallel do default(none)&
    !$omp private(rab, dr, vab, urab, rhoa, rhob, nwa, nwb, qa, qb, qc, Pa, Pb, n2wa, n2wb, j, i, r2) &
    !$omp private(projv, df, ddf, nlist, Jac, vba) &
    !$omp shared(dv, du, dh, dcf, n, pos, h, tt, v, den, c, p, om, mas, u, kcf, cf)&
    !$omp shared(dim, kr, kt, ptype, stepsize, gradv)
    do i = 1, n, stepsize
      ! print *, i, stepsize, ptype(i)
      ! read *
      if (ptype(i) /= 0) then
        dv(:,i) = 0.
        du(i) = 0.
        dh(i) = 0.
        dcf(i) = 0.
        call getneighbours(i, nlist)
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
          select case (tt)
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

            ! iderivtype = 1 (default): Brookshaw/Espanol-Revenga style derivatives
            !    del2v = \sum m_j/rho_j (v_i - v_j) * 2.0*abs(\nabla W_ij)/rij (h_i)
            ! iderivtype = 2: "standard" second derivatives with kernel
            !    del2v = \sum m_j/rho_j (v_j - v_i) nabla^2 W_ij (hi)
            vba(:) = v(:,j) - v(:,i)
            ! call get_kerntype(kt)
            ! if (kt == 1) then
              ! call get_jacobian(r, h(i), Jac)
              ! call get_n2w(r, h(i), n2wa)
              ! dv(:,i)  = dv(:,i) + mas(j)/den(j) * vba(:) * n2wa
              ! print*, vba(:)*n2wa
              ! dv(:,i)  = dv(:,i) + mas(j)/den(j) * vba(:)*(Jac(1,1) + Jac(2,2) + Jac(3,3))
              ! print*, vba(:)*(Jac(1,1) + Jac(2,2) + Jac(3,3))
              ! read*
            ! else
              call get_n2w(rab, h(i), n2wa)
              dv(:,i)  = dv(:,i) + mas(j)/den(j) * vba(:) * n2wa
            ! end if
          case(6)
            ! diff-graddiv
            ! iderivtype = 1 (default): Brookshaw/Espanol-Revenga style derivatives
            !    graddivv = \sum m_j/rho_j (5 vij.rij - vij) abs(\nabla W_ij)/rij (h_i)
            ! iderivtype = 2: "standard" second derivatives with kernel
            !    graddivv = \sum m_j/rho_j ((v_j - v_i).\nabla) \nabla W_ij
            call get_kerntype(kt)
            if (kt == 1) then
              ! call GradDivW(r, h(i), nwa)
              call get_jacobian(rab, h(i), Jac)
              dv(1,i) = dv(1,i) + mas(j)/den(j) * (vba(1)*Jac(1,1) + vba(2)*Jac(1,2) + vba(3)*Jac(1,3))
              dv(2,i) = dv(2,i) + mas(j)/den(j) * (vba(1)*Jac(2,1) + vba(2)*Jac(2,2) + vba(3)*Jac(2,3))
              dv(3,i) = dv(3,i) + mas(j)/den(j) * (vba(1)*Jac(3,1) + vba(2)*Jac(3,2) + vba(3)*Jac(3,3))
              ! ! dv(1,i)  = dv(1,i) + mas(j)/den(j) * vba(1)*Hes(1,1)
              ! ! dv(2,i)  = dv(2,i) + mas(j)/den(j) * vba(2)*Hes(2,2)
              ! ! dv(3,i)  = dv(3,i) + mas(j)/den(j) * vba(3)*Hes(3,3)
              ! print*, (vab(1)*Jac(2,1) + vab(2)*Jac(2,2) + vab(3)*Jac(2,3))
              !
              ! print*, vab(1)*(ddf*r(1)*r(2) - df*r(1)*r(2))/r2+&
              !         vab(2)*(ddf*r(2)*r(2)/r2 + df*(1*dr-r(2)*r(2))/r2)+ &
              !         vab(3)*(ddf*r(3)*r(2) - df*r(3)*r(2))/r2
              !
              ! print*, (ddf*r(2)*dot_product(vab(:),r(:)) + df*(vab(2)*dr - dot_product(vab(:),r(:))*r(2)))/r2
              !
              ! print*, "Diff: dvm - dvp: ", (vba(1)*Jac(2,1) + vba(2)*Jac(2,2) + vba(3)*Jac(2,3)) -&
              !                              (ddf*r(2)*projv/dr + df*(vab(2) - projv*r(2))/dr)
              !
              ! print*, "Diff: dvmx - dvp: ", vab(1)*(ddf*r(1)*r(2)-df*r(1)*r(2))/r2+&
              !         vab(2)*(ddf*r(2)*r(2)+df*(1*dr-r(2)*r(2)))/r2+ &
              !         vab(3)*(ddf*r(3)*r(2)-df*r(3)*r(2))/r2 -&
              !         (ddf*r(2)*dot_product(vab(:),r(:)) + df*(vab(2)*dr - dot_product(vab(:),r(:))*r(2)))/r2
              !

              ! read*
              ! print*, dv(:,i)
              ! dv(:,i) = 0.
              ! vab(:) = v(:,i) - v(:,j)
              ! projv = dot_product(vab(:),urab(:))
              ! call PureKernel(dr, h(i), df, ddf)
              ! ! graddivvi(:) = graddivvi(:) - pmass(j)*rho1j*(dr(:)*projv*grgrkerni + (dv(:) - projv*dr(:))*grkerni*rij1)
              ! dv(:,i) = dv(:,i) - mas(j)/den(j)*(r(:)/dr*projv*ddf + (vab(:) - projv*r(:)/dr)*df/dr)
              ! print*, dv(:,i)
              ! dv(:,i) = 0.
              ! dv(:,i) = dv(:,i) - mas(j)/den(j)*(r(:)/dr/dr*dot_product(vab(:),r(:))*ddf + (vab(:) - dot_product(vab(:),r(:))*r(:)/dr/dr)*df/h(i))
              ! print*, dv(:,i)
              ! dv(:,i) = 0.
              ! print*, '---------'
              ! read*
            elseif ( kt == 2 ) then
              ! Fab case
              ! urab(:) = r(:) / dr
              ! projv = dot_product(vab(:),urab(:))
              ! call get_n2w(r, h(i), n2wa)
              call get_jacobian(rab,h(i),Jac)
              ! dv(:,i) = dv(:,i) - 0.5*mas(j)/den(j) * ((dim + 2)*projv*urab(:) - vab(:)) * n2wa
              ! print*, dv(:,i)
              ! dv(:,i) = 0.
              dv(1,i) = dv(1,i) + mas(j)/den(j) * (vba(1)*Jac(1,1) + vba(2)*Jac(1,2) + vba(3)*Jac(1,3))
              dv(2,i) = dv(2,i) + mas(j)/den(j) * (vba(1)*Jac(2,1) + vba(2)*Jac(2,2) + vba(3)*Jac(2,3))
              dv(3,i) = dv(3,i) + mas(j)/den(j) * (vba(1)*Jac(3,1) + vba(2)*Jac(3,2) + vba(3)*Jac(3,3))
            elseif ( kt == 3 ) then
              call get_nw(rab, h(i), nwa)
              dv(1,i) = dv(1,i) + mas(j)/den(j) * 0.5*((2*gradv(1,1,j)-2*gradv(1,1,i))*nwa(1) + &
                            (gradv(1,2,j) + gradv(2,1,j) - gradv(1,2,i) - gradv(2,1,i))*nwa(2) + &
                            (gradv(1,3,j) + gradv(3,1,j) - gradv(1,3,i) - gradv(3,1,i))*nwa(3))
              dv(2,i) = dv(2,i) + mas(j)/den(j) * 0.5*((gradv(2,1,j) + gradv(1,2,j) - gradv(2,1,i) - gradv(1,2,i))*nwa(1) + &
                            (2*gradv(2,2,j) - 2*gradv(2,2,i))*nwa(2) + &
                            (gradv(2,3,j) + gradv(3,2,j) - gradv(2,3,i) - gradv(3,2,i))*nwa(3))
              dv(3,i) = dv(3,i) + mas(j)/den(j) * 0.5*((gradv(3,1,j) + gradv(1,3,j) - gradv(3,1,i) - gradv(1,3,i))*nwa(1) + &
                            (gradv(3,2,j) + gradv(2,3,j) - gradv(3,2,i) - gradv(2,3,i))*nwa(2) + &
                            (2*gradv(3,3,j) - 2*gradv(3,3,i))*nwa(3))
            end if
          case default
            print *, 'Task type was not defined in circuit2 inside circle'
            stop
          end select
          if ( i == 88001) then
            ! print*, i, j, dv(:,i)
            ! read*
          end if
        end do
        select case (tt)
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

  subroutine fungradient(dim, t, m, d, x, v, h, nv)
    real, allocatable, intent(in)  :: m(:), d(:), x(:,:), v(:,:), h(:)
    integer, allocatable, intent(in) :: t(:)
    real, allocatable, intent(out) :: nv(:,:,:)
    integer, intent(in)            :: dim

    integer, allocatable :: nlist(:)
    integer :: n, i, j, l, ni, nj, li
    real    :: vba(3), nw(3), rab(3)

    n = size(m)
    allocate(nv(3,3,n))

    !$omp parallel do default(none)&
    !$omp private(rab, vba, nw, i, j, l, ni, nj, li, nlist) &
    !$omp shared(dim, t, m, d, x, v, h, nv, n, stepsize)
    do i = 1,n
      nv(:,:,i) = 0.
        call getneighbours(i, nlist)
        do l = 1, size(nlist)
        j = nlist(l)
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
  end subroutine fungradient
end module circuit2
