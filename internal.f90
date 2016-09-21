module internal
  use kernel
  use printer
  use setup

 implicit none

 public :: get_kinetic_energy, derivs

 private

contains

  subroutine get_density(n, pos, mas, den, slen)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n), mas(n), slen(n)
    real, intent(out)   :: den(n)
    real                :: w, dw, q
    integer             :: i, j

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

  subroutine equation_of_state(n, den, u, P, c, gamma)
    integer, intent(in) :: n
    real, intent(in)    :: den(n), u(n), c, gamma
    real, intent(out)   :: P(n)
    integer             :: i

    do i = 1, n
      P(i) = den(i) * c * c
      ! P(i) = (gamma - 1) * den(i) * u(i)
    end do
  end subroutine equation_of_state

  subroutine get_accel(n, c, pos, vel, mas, den, slen, P, acc, du)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n), mas(n), slen(n), den(n), P(n), vel(n), c
    real, intent(out)   :: acc(n), du(n)
    real                :: Pi, Pj, wi, wj, dwi, dwj, qi, qj, dx, qa, qb
    integer             :: i, j

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
            call art_viscosity(den(i), den(j), vel(i), vel(j), pos(i), pos(j), c, qa, qb)
            ! v + 2 in dimmentions
            dwi = dwi * (pos(i) - pos(j)) / slen(i) ** 3
            dwj = dwj * (pos(i) - pos(j)) / slen(j) ** 3
            Pi = (P(i) + qa) * dwi / den(i)**2
            Pj = (P(j) + qb) * dwj / den(j)**2
            acc(i) = acc(i) - mas(j) * (Pi + Pj)
            du(i) = du(i) + mas(j) * (vel(i) - vel(j)) * Pi
          end if
        end if
      end do
    end do
  end subroutine get_accel

  subroutine art_viscosity(da, db, va, vb, ra, rb, c, qa, qb)
    real, intent(in)  :: da, db, va, vb, ra, rb, c
    real, intent(out) :: qa, qb
    real              :: alpha, betta, vab, vba, r
    qa = 0.
    qb = 0.
    alpha = 1.
    betta = 2.

    vab = va - vb
    vba = -vab
    r = abs(ra - rb)
    if (vab * r < 0) then
      qa = -0.5 * da * (alpha*c - betta*(vab * r)) * (vab * r)
    end if
    if (vba * r < 0) then
      qb = -0.5 * db * (alpha*c - betta*(vba * r)) * (vba * r)
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

  subroutine derivs(t, n, bn, pos, vel, mas, den, slen, pres, acc, u, du, sos, sk)
    integer, intent(in)           :: n, bn
    real, intent(in)              :: sos, sk
    real, intent(out)             :: pos(n), vel(n), mas(n), den(n), slen(n), pres(n), acc(n), u(n), du(n)
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

    call equation_of_state(n, den, u, pres, sos, 1.2)

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
