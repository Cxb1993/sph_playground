module internal
  use kernel
  use printer
  use setup
  use eos
  use omp_lib

 implicit none

 public :: derivs

 private

contains

  subroutine get_density(n, pos, mas, sk, sln, den, om)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n,3), mas(n), sk
    real, intent(out)   :: sln(n), den(n), om(n)
    real                :: w, dwdh, r(3), dr, dnrfh, nrfh, h0, hp, hn
    integer             :: i, j, dim, itern

    call get_dim(dim)

    !$OMP PARALLEL
    !$OMP DO PRIVATE(r, dr, dwdh, w, dnrfh, nrfh, h0, hp, hn, itern)
    do i = 1, n
      den(i) = 0.
      om(i) = 0.
      j = 0
      if (i.ne.j) then
        do j = 1, n
          r(:) = pos(i,:) - pos(j,:)
          dr = sqrt(dot_product(r(:),r(:)))
          if (dr < 2. * sln(i)) then
            call get_dw_dh(dr, sln(i), dwdh)
            call get_w(dr, sln(i), w)
            den(i) = den(i) + mas(j) * w
            om(i) = om(i) + mas(j) * dwdh
          end if
        end do
      end if
      om(i) = 1. - om(i) * (- sln(i) / (3. * den(i)))
      h0 = sln(i)
      hp = 1000 * h0
      hn = h0
      itern = 50
      do while ((abs(hn - hp) / h0 > 0.0001).and.(itern.gt.0))
        hp = hn
        dnrfh = - 3. * den(i) * om(i) / hp
        nrfh  = mas(i) * (sk / hp) ** dim - den(i)
        hn = hp - nrfh / dnrfh
        itern = itern - 1
        if (itern.eq.0) then
          print *, 'Something wrong it softening length. i:', i, 'h_new', hn
          print *, pos(i,:)
          print *, hp, nrfh, dnrfh
          stop
        end if
      end do
      sln(i) = hn
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine get_density

  subroutine get_accel(n, c, pos, vel, acc, mas, den, sln, om, P, u, du, dh)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n,3), vel(n,3), mas(n), sln(n), den(n), P(n), c(n), om(n)
    real, intent(out)   :: acc(n,3), u(n), du(n), dh(n)
    real                :: dr, di, dj, qa, qb, qc
    real                :: nwi(3), nwj(3), r(3), vab(3), urab(3), Pi(3), Pj(3)
    integer             :: i, j, dim

    call get_dim(dim)
    !$OMP PARALLEL
    !$OMP DO PRIVATE(r, dr, vab, urab,di, dj, nwi, nwj, qa, qb, qc, Pi, Pj)
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
            call art_viscosity(di, dj, vab, urab, c(i), c(j), qa, qb)
            call art_termcond(P(i), P(j), di, dj, qc)
            Pi(:) = (P(i) + qa) * nwi(:) / (di**2 * om(i))
            Pj(:) = (P(j) + qb) * nwj(:) / (dj**2 * om(j))

            acc(i,:) = acc(i,:) - mas(j) * (Pi(:) + Pj(:))

            du(i) = du(i) + mas(j) * dot_product(vab(:),Pi(:)) &
                          + mas(j) / (0.5 *(di + dj)) * qc * (u(i) - u(j)) * &
                          0.5 * dot_product((nwi(:) + nwj(:)),urab(:))

            dh(i) = dh(i) + mas(j) * dot_product(vab(:), nwi(:))
          end if
        end if
      end do
      dh(i) =  (- sln(i) / (dim * den(i))) * dh(i) / om(i)
    end do
    !$OMP END DO
    !$OMP END PARALLEL

  end subroutine get_accel

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

  subroutine derivs(t, n, pos, vel, acc, mas, den, h, dh, om, prs, c, uei, due, sos, sk, gamma)
    integer, intent(in)           :: n
    real, intent(in)              :: sos, sk, gamma
    real, intent(inout)           :: pos(n,3), vel(n,3), acc(n,3)
    real, intent(inout)           :: mas(n), den(n), h(n), dh(n), prs(n), c(n), uei(n), due(n), om(n)
    character (len=*), intent(in) :: t

    call get_density(n, pos, mas, sk, h, den, om)
    if (t.eq.'periodic') then
      call set_periodic(n, den)
    end if

    select case (t)
      case ('periodic')
        call eos_isothermal(n, den, prs, sos)
      case ('shock_fixed')
        call eos_adiabatic(n, den, uei, prs, c, gamma)
    end select

    if (t.eq.'periodic') then
      call set_periodic(n, prs)
    end if

    call get_accel(n, c, pos, vel, acc, mas, den, h, om, prs, uei, due, dh)

    select case (t)
      case ('periodic')
        call set_periodic(n, acc)
      case ('shock_fixed')
        call set_fixed3(acc)
        call set_fixed1(due)
      end select

  end subroutine derivs

end module internal
