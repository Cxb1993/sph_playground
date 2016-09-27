module internal
  use kernel
  use printer
  use setup

 implicit none

 public :: get_kinetic_energy, derivs, get_denergy

 private

contains

  subroutine get_density(n, pos, mas, den, slen)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n), mas(n), slen(n)
    real, intent(out)   :: den(n)
    real                :: w, dw, q
    integer             :: i, j
    j = 0

    do i = 1, n
      den(i) = 0.
      if (i /= j) then
        do j = 1, n
           q = abs(pos(i)-pos(j))/slen(i)
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

  subroutine get_accel(n, c, pos, vel, mas, den, slen, P, acc, du)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n), mas(n), slen(n), den(n), P(n), vel(n), c(n)
    real, intent(out)   :: acc(n), du(n)
    real                :: Pi, Pj, wi, wj, dwi, dwj, qi, qj, dx, qa, qb
    integer             :: i, j

    qa = 0.
    qb = 0.

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
            ! v + 2 in dimmentions
            dwi = dwi * (pos(i) - pos(j)) / slen(i) ** 3
            dwj = dwj * (pos(i) - pos(j)) / slen(j) ** 3
            Pi = (P(i) + qa) * dwi / den(i)**2
            Pj = (P(j) + qb) * dwj / den(j)**2
            acc(i) = acc(i) - mas(j) * (Pi + Pj)
            du(i) = du(i) + mas(j) * (vel(i) - vel(j)) * Pi !&
                          ! + mas(j) * vsigu * (en(i)-en(j)) * 0.5(dwi/den(i) + dwj/den(j))
          end if
        end if
      end do
    end do
  end subroutine get_accel

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
      ! print *, qa, qb
      ! read *
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

  subroutine get_denergy(n, mas, vel, acc, du, de)
    integer, intent(in) :: n
    real, intent(in)    :: mas(n), vel(n), acc(n), du(n)
    real, intent(out)   :: de
    integer             :: i

    de = 0.
    do i = 1, n
      de = de + mas(i) * ( vel(i) * acc(i) + du(i))
    end do
  end subroutine get_denergy

  subroutine get_slength(n, mas, den, slen, sk)
    integer, intent(in) :: n
    real, intent(in)    :: sk, mas(n), den(n)
    real, intent(out)   :: slen(n)
    integer             :: i

    do i = 1, n
      slen(i) = sk * (mas(i) / den(i))
    end do
  end subroutine get_slength

  subroutine derivs(t, n, bn, pos, vel, mas, den, slen, pres, c, acc, u, du, sos, sk, gamma)
    integer, intent(in)           :: n, bn
    real, intent(in)              :: sos, sk, gamma
    real, intent(out)             :: pos(n), vel(n), mas(n), den(n), slen(n), pres(n), c(n), acc(n), u(n), du(n)
    real                          :: de
    character (len=*), intent(in) :: t
    integer                       :: i

    do i = 1, 3
      call get_density(n, pos, mas, den, slen)

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

    call get_accel(n, c, pos, vel, mas, den, slen, pres, acc, du)

    select case (t)
      case ('periodic')
        call set_periodic(n, bn, acc)
      case ('fixed')
        call set_fixed(n, bn, acc)
        call set_fixed(n, bn, du)
    end select
    ! call get_denergy(n, mas, vel, acc, du, de)
    ! print *, de
    ! read *

  end subroutine derivs

end module internal
