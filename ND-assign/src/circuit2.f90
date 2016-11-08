module circuit2
  use omp_lib
  use kernel
  use BC

  implicit none

  public :: make_c2

  private

contains

  subroutine make_c2(n, c, pos, vel, acc, mas, den, sln, om, P, u, du, dh, cf, dcf, kcf)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n,3), vel(n,3), mas(n), sln(n), den(n), P(n), c(n), om(n),&
                           cf(n), kcf(n)
    real, intent(out)   :: acc(n,3), u(n), du(n), dh(n), dcf(n)
    real                :: dr, di, dj, qa, qb, qc, n2w, fab
    real                :: nwi(3), nwj(3), r(3), vab(3), urab(3), Pi(3), Pj(3)
    integer             :: i, j, dim
    character(len=40)   :: t

    call get_dim(dim)
    call get_tasktype(t)

    !$OMP PARALLEL
    !$OMP DO PRIVATE(r, dr, vab, urab, di, dj, nwi, nwj, qa, qb, qc, Pi, Pj, n2w, fab)
    do i = 1, n
      acc(i,:) = 0.
      du(i) = 0.
      dh(i) = 0.
      dcf(i) = 0.
      do j = 1, n
        if (i.ne.j) then
          r(:) = pos(i,:) - pos(j,:)
          dr = sqrt(dot_product(r(:),r(:)))
          if ((dr < 3. * sln(i)) .or. (dr < 3. * sln(j))) then
            select case (t)
            case ('hydroshock')
              qa = 0.
              qb = 0.
              qc = 0.
              vab(:) = vel(i,:) - vel(j,:)
              urab(:) = r(:) / dr
              di = den(i)
              dj = den(j)

              call get_nw(r, sln(i), nwi)
              call get_nw(r, sln(j), nwj)
              call art_viscosity(di, dj, vab, urab, c(i), c(j), qa, qb)
              call art_termcond(P(i), P(j), di, dj, qc)
              Pi(:) = (P(i) + qa) * nwi(:) / (di**2 * om(i))
              Pj(:) = (P(j) + qb) * nwj(:) / (dj**2 * om(j))

              acc(i,:) = acc(i,:) - mas(j) * (Pi(:) + Pj(:))

              du(i)   = du(i) + mas(j) * dot_product(vab(:),Pi(:)) &
                            + mas(j) / (0.5 *(di + dj)) * qc * (u(i) - u(j)) * &
                            0.5 * dot_product((nwi(:) + nwj(:)),urab(:))

              dh(i)   = dh(i) + mas(j) * dot_product(vab(:), nwi(:))
            case ('heatslab')
              call get_n2w(dr, sln(i), n2w)
              call get_Fab(r, sln(i), fab)
              ! 2 * kcf(i) * kcf(j) /
              du(i) = du(i) - mas(j) / (den(i) * den(j)) * 2 * kcf(i) * kcf(j) / (kcf(i) + kcf(j)) &
                            * (cf(i) - cf(j)) * n2w

              ! print *, (-0.13 + 2.22/(50*(dr/sln(i) - 1.12)**2 + 1)), &
              !          Fab, &
              !          n2w
              ! read *
            end select
          end if
        end if
      end do
      dh(i) =  (- sln(i) / (dim * den(i))) * dh(i) / om(i)
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine make_c2

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
