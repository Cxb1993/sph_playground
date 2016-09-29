module internal
  use kernel
  use printer
  use setup

 implicit none

 public :: get_kinetic_energy, derivs

 private

contains

  subroutine get_density(n, pos, mas, den, sln)
    integer, intent(in) :: n
    real, intent(in)    :: pos(3,n), mas(n), sln(n)
    real, intent(out)   :: den(n)
    real                :: w, dw, q
    integer             :: i, j
    j = 0

    do i = 1, n
      den(i) = 0.
      if (i /= j) then
        do j = 1, n
           q = sqrt((pos(1,i)-pos(1,j))**2 + (pos(2,i)-pos(2,j))**2 + (pos(3,i)-pos(3,j))**2 )/slen(i)
           if (q < 2.) then
              call get_kernel(q, w, dw)
              den(i) = den(i) + mas(j) * w / slen(i)
           endif
        end do
      end if
    end do
  end subroutine get_density

  subroutine eos_adiabatic(n, den, u, P, c, gamma)
    integer, intent(in) :: n
    real, intent(in)    :: den(n), u(n), gamma
    real, intent(out)   :: P(n), c(n)
    integer             :: i

    do i = 1, n
      P(i) = (gamma - 1) * den(i) * u(i)
      c(i) = sqrt(gamma * P(i) / den(i))
    end do
  end subroutine eos_adiabatic

  subroutine eos_isothermal(n, den, P, c)
    integer, intent(in) :: n
    real, intent(in)    :: den(n), c
    real, intent(out)   :: P(n)
    integer             :: i

    do i = 1, n
      P(i) = den(i) * c * c
    end do
  end subroutine eos_isothermal

  subroutine get_accel(n, c, pos, vel, mas, den, slen, P, acc, u, du)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n), mas(n), slen(n), den(n), P(n), vel(n), c(n)
    real, intent(out)   :: acc(n), u(n), du(n)
    real                :: Pi, Pj, wi, wj, dwi, dwj, qi, qj, dx, qa, qb, qc
    integer             :: i, j

    qa = 0.
    qb = 0.
    qc = 0.

    do i = 1, n
      acc(i) = 0.
      du(i) = 0.
      do j = 1, n
        if (i /= j) then
          dx = abs(pos(i)-pos(j))
          qi = dx/slen(i)
          qj = dx/slen(j)
          if (qi < 2. .or. qj < 2.) then
            call get_kernel(qi, wi, dwi)
            call get_kernel(qj, wj, dwj)
            call art_viscosity(den(i), den(j), vel(i), vel(j), pos(i), pos(j), c(i), c(j), qa, qb)
            call art_termcond(P(i), P(j), den(i), den(j), pos(i), pos(j), qc)
            ! v + 2 in dimmentions
            dwi = dwi * (pos(i) - pos(j)) / slen(i) ** 3
            dwj = dwj * (pos(i) - pos(j)) / slen(j) ** 3
            Pi = (P(i) + qa) * dwi / den(i)**2
            Pj = (P(j) + qb) * dwj / den(j)**2
            acc(i) = acc(i) - mas(j) * (Pi + Pj)
            du(i) = du(i) + mas(j) * (vel(i) - vel(j)) * Pi &
                          + mas(j) / (0.5 *(den(i) + den(j))) * qc * (u(i) - u(j)) * &
                          0.5 * (dwi + dwj) * ((pos(i) - pos(j))/abs(pos(i) - pos(j)))
          end if
        end if
      end do
    end do
  end subroutine get_accel

  subroutine art_termcond(pa, pb, da, db, ra, rb, vsigu)
    real, intent(in)  :: pa, pb, da, db, ra, rb
    real, intent(out) :: vsigu

    vsigu = sqrt(abs(pa - pb)/(0.5 * (da + db)))

  end subroutine art_termcond

  subroutine art_viscosity(da, db, va, vb, ra, rb, c, qa, qb)
    real, intent(in)  :: da, db, va, vb, ra, rb, c
    real, intent(out) :: qa, qb
    real              :: alpha, betta, vab, vba, rab, rba
    qa = 0.
    qb = 0.
    alpha = 1.
    betta = 2.

    vab = va - vb
    vba = -vab
    rab = ra - rb
    rba = -rab
    if (vab * rab < 0) then
      qa = -0.5 * da * (alpha*c - betta*(vab * rab)) * (vab * rab)
    end if
    if (vba * rba < 0) then
      qb = -0.5 * db * (alpha*c - betta*(vba * rba)) * (vba * rba)
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

  subroutine get_slength(n, mas, den, slen, sk)
    integer, intent(in) :: n
    real, intent(in)    :: sk, mas(n), den(n)
    real, intent(out)   :: slen(n)
    integer             :: i

    do i = 1, n
      slen(i) = sk * (mas(i) / den(i))
    end do
  end subroutine get_slength

  subroutine derivs(t, n, bn, pos, vel, acc, mas, den, sln, prs, c, uei, due, sos, sk, gamma)
    integer, intent(in)           :: n, bn
    real, intent(in)              :: sos, sk, gamma
    real, intent(out)             :: pos(3,n), vel(n), acc(n), mas(n), den(n), sln(n), prs(n), c(n), uei(n), due(n)
    character (len=*), intent(in) :: t
    integer                       :: i

    do i = 1, 3
      call get_density(n, pos, mas, den, sln)

      if (t .EQ. 'periodic') then
        call set_periodic(n, bn, den)
      end if

      call get_slength(n, mas, den, slen, sk)

      if (t .EQ. 'periodic') then
        call set_periodic(n, bn, slen)
      end if
    end do

    select case (t)
      case ('periodic')
        call eos_isothermal(n, den, pres, sos)
      case ('fixed')
        call eos_adiabatic(n, den, u, pres, c, gamma)
    end select

    if (t .EQ. 'periodic') then
      call set_periodic(n, bn, pres)
    end if

    call get_accel(n, sos, pos, vel, mas, den, slen, pres, acc, du)

    select case (t)
      case ('periodic')
        call set_periodic(n, bn, acc)
      case ('fixed')
        call set_fixed(n, bn, vel)
    end select
  end subroutine derivs

end module internal
