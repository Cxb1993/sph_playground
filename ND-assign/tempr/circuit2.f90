module tempr_circuit2
  use omp_lib
  use kernel

  implicit none

  public :: tempr_hydro_c2, tempr_solid_c2

  private

contains

  subroutine tempr_hydro_c2(n, c, pos, vel, acc, mas, den, sln, om, P, u, du, dh, cf, dcf, kcf)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n,3), vel(n,3), mas(n), sln(n), den(n), P(n), c(n), om(n), cf(n)
    real, intent(out)   :: acc(n,3), u(n), du(n), dh(n), dcf(n), kcf(n)
    real                :: dr, di, dj, qa, qb, qc, n2y
    real                :: nwi(3), nwj(3), r(3), vab(3), urab(3), Pi(3), Pj(3)
    integer             :: i, j, dim

    call get_dim(dim)
    ! print *, 'Dim in circuit2: ', dim
    ! read *
    !$OMP PARALLEL
    !$OMP DO PRIVATE(r, dr, vab, urab,di, dj, nwi, nwj, qa, qb, qc, Pi, Pj, n2y)
    do i = 1, n
      acc(i,:) = 0.
      du(i) = 0.
      dh(i) = 0.
      do j = 1, n
        if (i.ne.j) then
          r(:) = pos(i,:) - pos(j,:)
          dr = sqrt(dot_product(r(:),r(:)))
          if (dr < 2. * sln(i) .or. dr < 2. * sln(j)) then
            qa = 0.
            qb = 0.
            qc = 0.
            vab(:) = vel(i,:) - vel(j,:)
            urab(:) = r(:) / dr
            di = den(i)
            dj = den(j)

            call get_nabla_w(r, sln(i), nwi)
            call get_nabla_w(r, sln(j), nwj)
            call get_n2y(dr, sln(i), n2y)
            call art_viscosity(di, dj, vab, urab, c(i), c(j), qa, qb)
            call art_termcond(P(i), P(j), di, dj, qc)
            Pi(:) = (P(i) + qa) * nwi(:) / (di**2 * om(i))
            Pj(:) = (P(j) + qb) * nwj(:) / (dj**2 * om(j))

            acc(i,:) = acc(i,:) - mas(j) * (Pi(:) + Pj(:))

            du(i)   = du(i) + mas(j) * dot_product(vab(:),Pi(:)) &
                          + mas(j) / (0.5 *(di + dj)) * qc * (u(i) - u(j)) * &
                          0.5 * dot_product((nwi(:) + nwj(:)),urab(:))

            dh(i)   = dh(i) + mas(j) * dot_product(vab(:), nwi(:))

            dcf(i) = du(i) - 4 * mas(j) / (di * dj) * kcf(i) * kcf(j) / (kcf(i) + kcf(j)) &
                             * (cf(i) - cf(j)) * n2y * 0.01 ** 2
          end if
        end if
      end do
      dh(i) =  (- sln(i) / (dim * den(i))) * dh(i) / om(i)
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine tempr_hydro_c2

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

  subroutine tempr_solid_c2(n, pos, mas, den, sln, u, du, kcf)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n,3), mas(n), sln(n), den(n), u(n), kcf(n)
    real, intent(out)   :: du(n)
    real                :: dr, n2y, r(3), nw(3), Fab
    integer             :: i, j

    !$OMP PARALLEL
    !$OMP DO PRIVATE(r, dr, n2y, nw, Fab)
    do i = 1, n
      du(i) = 0.
      do j = 1, n
        if (i.ne.j) then
          r(:) = pos(i,:) - pos(j,:)
          dr = sqrt(dot_product(r(:),r(:)))
          if (dr < 2. * sln(i)) then
            call get_n2y(dr, sln(i), n2y)
            call get_nabla_w(r, sln(i), nw)
            Fab = sqrt(dot_product(nw,nw)) / dr
            ! print *, sqrt(dot_product(nw,nw)) / dr, n2y
            ! read *

            du(i) = du(i) - 4. * mas(j) / (den(i) * den(j)) * kcf(i) * kcf(j) / (kcf(i) + kcf(j)) &
                          * (u(i) - u(j)) * Fab
          end if
        end if
      end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine tempr_solid_c2
end module tempr_circuit2
