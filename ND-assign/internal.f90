module internal
  use kernel
  use printer
  use setup
  use eos

 implicit none

 public :: get_kinetic_energy, derivs

 private

contains

  subroutine get_density(n, pos, mas, sln, den, om)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n,3), mas(n), sln(n)
    real, intent(out)   :: den(n), om(n)
    real                :: w, dwdh, r(3), dr
    integer             :: i, j
    j = 0

    do i = 1, n
      den(i) = 0.
      om(i) = 0.
      if (i.ne.j) then
        do j = 1, n
          r(:) = pos(i,:) - pos(j,:)
          dr = sqrt(dot_product(r(:),r(:)))
          if (dr < 2. * sln(i)) then
            call get_dw_dh(dr, sln(i), w, dwdh)
            den(i) = den(i) + mas(j) * w
            om(i) = om(i) + mas(j) * dwdh
          end if
        end do
      end if
      om(i) = 1. - om(i) * (- sln(i) / (3. * den(i)))
    end do
  end subroutine get_density

  subroutine get_accel(n, c, pos, vel, acc, mas, den, sln, om, P, u, du)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n,3), mas(n), sln(n), den(n), P(n), vel(n), c(n), om(n)
    real, intent(out)   :: acc(n), u(n), du(n)
    real                :: Pi, Pj, wi, wj, nwi(3), nwj(3), dr, qa, qb, qc, r(3)
    integer             :: i, j

    qa = 0.
    qb = 0.
    qc = 0.

    do i = 1, n
      acc(i) = 0.
      du(i) = 0.
      do j = 1, n
        if (i.ne.j) then
          r(:) = pos(i,:) - pos(j,:)
          dr = sqrt(dot_product(r(:),r(:)))
          if (dr < 2. * sln(i) .or. dr < 2. * sln(j)) then
            call get_nabla_w(r, sln(i), wi, nwi)
            call get_nabla_w(r, sln(j), wj, nwj)
            call art_viscosity(den(i), den(j), vel(i), vel(j), pos(i,1), pos(j,1), c(i), c(j), qa, qb)
            call art_termcond(P(i), P(j), den(i), den(j), qc)
            ! v + 2 in dimmentions
            Pi = (P(i) + qa) * nwi(1) / (den(i)**2 * om(i))
            Pj = (P(j) + qb) * nwj(1) / (den(j)**2 * om(j))
            acc(i) = acc(i) - mas(j) * (Pi + Pj)
            du(i) = du(i) + mas(j) * (vel(i) - vel(j)) * Pi &
                          + mas(j) / (0.5 *(den(i) + den(j))) * qc * (u(i) - u(j)) * &
                          0.5 * (nwi(1) + nwj(1)) * ((pos(i,1) - pos(j,1))/abs(pos(i,1) - pos(j,1)))
          end if
        end if
      end do
    end do
  end subroutine get_accel

  subroutine art_termcond(pa, pb, da, db, vsigu)
    real, intent(in)  :: pa, pb, da, db
    real, intent(out) :: vsigu

    vsigu = sqrt(abs(pa - pb)/(0.5 * (da + db)))
  end subroutine art_termcond

  subroutine art_viscosity(da, db, va, vb, ra, rb, ca, cb, qa, qb)
    real, intent(in)  :: da, db, va, vb, ra, rb, ca, cb
    real, intent(out) :: qa, qb
    real              :: alpha, betta, vab, rab
    qa = 0.
    qb = 0.
    alpha = 1.
    betta = 2.

    vab = va - vb
    rab = (ra - rb)/abs(ra - rb)
    if (vab * rab < 0) then
      qa = -0.5 * da * (alpha*ca - betta*(vab * rab)) * (vab * rab)
      qb = -0.5 * db * (alpha*cb - betta*(vab * rab)) * (vab * rab)
    end if
  end subroutine art_viscosity

  subroutine get_kinetic_energy(n, mas, vel, e)
    integer, intent(in) :: n
    real, intent(in)    :: mas(n), vel(n)
    real, intent(out)   :: e
    integer             :: i

    e = 0.
    do i = 1, n
      e = e + 0.5 * mas(i) * vel(i) ** 2
    end do
  end subroutine get_kinetic_energy

  subroutine get_slength(n, mas, den, sln, sk)
    integer, intent(in) :: n
    real, intent(in)    :: sk, mas(n), den(n)
    real, intent(out)   :: sln(n)
    integer             :: i

    do i = 1, n
      sln(i) = sk * (mas(i) / den(i))
    end do
  end subroutine get_slength

  subroutine derivs(t, n, bn, pos, vel, acc, mas, den, sln, om, prs, c, uei, due, sos, sk, gamma)
    integer, intent(in)           :: n, bn
    real, intent(in)              :: sos, sk, gamma
    real, intent(out)             :: pos(n,3), vel(n), acc(n), mas(n), den(n), sln(n), prs(n), c(n), uei(n), due(n), om(n)
    character (len=*), intent(in) :: t
    integer                       :: i

    do i = 1, 3
      call get_density(n, pos, mas, sln, den, om)

      if (t.eq.'periodic') then
        call set_periodic(n, bn, den)
      end if

      call get_slength(n, mas, den, sln, sk)

      if (t.eq.'periodic') then
        call set_periodic(n, bn, sln)
      end if
    end do

    select case (t)
      case ('periodic')
        call eos_isothermal(n, den, prs, sos)
      case ('shock_fixed')
        call eos_adiabatic(n, den, uei, prs, c, gamma)
    end select

    if (t.eq.'periodic') then
      call set_periodic(n, bn, prs)
    end if

    call get_accel(n, c, pos, vel, acc, mas, den, sln, om, prs, uei, due)

    select case (t)
      case ('periodic')
        call set_periodic(n, bn, acc)
      case ('shock_fixed')
        call set_fixed(n, bn, acc)
        call set_fixed(n, bn, due)
      end select
  end subroutine derivs

end module internal
