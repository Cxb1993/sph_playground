module circuit2
  use omp_lib
  use kernel
  use BC

  implicit none

  public :: c2

  private

contains

  subroutine c2(n, c, pos, vel, acc, mas, den, sln, om, P, u, du, dh, cf, dcf, kcf)
    integer, intent(in) :: n
    real, intent(in)    :: pos(3,n), vel(3,n), mas(n), sln(n), den(n), P(n), c(n), om(n),&
                           cf(n), kcf(n)
    real, intent(out)   :: acc(3,n), u(n), du(n), dh(n), dcf(n)
    real                :: dr, di, dj, qa, qb, qc, n2w, kr, r2
    real                :: nwi(3), nwj(3), r(3), vab(3), urab(3), Pi(3), Pj(3)
    integer             :: i, j, dim
    integer             :: tt!, kt

    call get_dim(dim)
    call get_tasktype(tt)
    call get_krad(kr)

    !$omp parallel do default(none)&
    !$omp private(r, dr, vab, urab, di, dj, nwi, nwj, qa, qb, qc, Pi, Pj, n2w, j, i, r2) &
    !$omp shared(acc, du, dh, dcf, n, pos, sln, tt, vel, den, c, p, om, mas, u, kcf, cf, dim, kr)
    do i = 1, n
      acc(:,i) = 0.
      du(i) = 0.
      dh(i) = 0.
      dcf(i) = 0.
      do j = 1, n
        if (i.ne.j) then
          r(:) = pos(:,i) - pos(:,j)
          r2 = dot_product(r(:),r(:))
          if ((r2 < (kr * sln(i))**2).or.(r2 < (kr * sln(j))**2)) then
          ! if (r2 < (kr * sln(i))**2) then
            dr = sqrt(r2)
            select case (tt)
            case (1)
              qa = 0.
              qb = 0.
              qc = 0.
              vab(:) = vel(:,i) - vel(:,j)
              urab(:) = r(:) / dr
              di = den(i)
              dj = den(j)

              call get_nw(r, sln(i), nwi)
              call get_nw(r, sln(j), nwj)
              call get_n2w(r, sln(i), n2w)

              call art_viscosity(di, dj, vab, urab, c(i), c(j), qa, qb)
              call art_termcond(P(i), P(j), di, dj, qc)
              Pi(:) = (P(i) + qa) * nwi(:) / (di**2 * om(i))
              Pj(:) = (P(j) + qb) * nwj(:) / (dj**2 * om(j))

              acc(:,i) = acc(:,i) - mas(j) * (Pi(:) + Pj(:))

              du(i)   = du(i) + mas(j) * dot_product(vab(:),Pi(:)) &
                            + mas(j) / (0.5 *(di + dj)) * qc * (u(i) - u(j)) * &
                            0.5 * dot_product((nwi(:) + nwj(:)),urab(:))

              dh(i)   = dh(i) + mas(j) * dot_product(vab(:), nwi(:))
            case (2, 3)
              call get_n2w(r, sln(i), n2w)

              du(i) = du(i) - mas(j) / (den(i) * den(j)) * 2. * kcf(i) * kcf(j) &
                      / (kcf(i) + kcf(j)) * (cf(i) - cf(j)) * n2w
              ! du(i) = du(i) - mas(j) / (den(i) * den(j)) * (kcf(i) + kcf(j)) / 2. &
              !               * (cf(i) - cf(j)) * n2w
            end select
          end if
        end if
      end do
      dh(i) =  (- sln(i) / (dim * den(i))) * dh(i) / om(i)
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
end module circuit2
