module internal
  use kernel
  use printer
  use setup

 implicit none

 public :: get_kinetic_energy, derivs

 private

contains

  subroutine get_density(n, bn, pos, mas, den, slen)
    integer, intent(in) :: n, bn
    real, intent(in)    :: pos(n), mas(n), slen(n)
    real, intent(out)   :: den(n)
    real                :: w, dw, q
    integer             :: i, j

    do i = bn, n - bn
      den(i) = 0.
      do j = 1, n
         q = abs(pos(i)-pos(j))/slen(i)
         if (q < 2.) then
            call get_kernel(q, w, dw)
            den(i) = den(i) + mas(j) * w / slen(i)
         endif
      end do
    end do
    call set_periodic(n, bn, den)
  end subroutine get_density

  subroutine equation_of_state(n, bn, den, P, c)
    integer, intent(in) :: n, bn
    real, intent(in)    :: den(n), c
    real, intent(out)   :: P(n)
    integer             :: i

    do i = bn, n - bn
      P(i) = den(i) * c * c
    end do
    call set_periodic(n, bn, P)
  end subroutine equation_of_state

  subroutine get_accel(n, bn, pos, mas, den, slen, P, acc)
    integer, intent(in) :: n, bn
    real, intent(in)    :: pos(n), mas(n), slen(n), den(n), P(n)
    real, intent(out)   :: acc(n)
    real                :: Pi, Pj, wi, wj, dwi, dwj, qi, qj, dx
    integer             :: i, j

    do i = bn, n - bn
      acc(i) = 0.
      do j = 1, n
        if (i /= j) then
          dx = abs(pos(i)-pos(j))
          qi = dx/slen(i)
          qj = dx/slen(j)
          if (qi < 2. .or. qj < 2.) then
            call get_kernel(qi, wi, dwi)
            call get_kernel(qj, wj, dwj)
            dwi = dwi * (pos(i) - pos(j)) / slen(i) ** 2
            dwj = dwj * (pos(i) - pos(j)) / slen(j) ** 2
            Pi = P(i) * dwi / (den(i) * den(i))
            Pj = P(j) * dwj / (den(j) * den(j))
            acc(i) = acc(i) - mas(j) * (Pi + Pj)
          end if
        end if
      end do
    end do
    call set_periodic(n, bn, acc)
  end subroutine get_accel

  subroutine get_kinetic_energy(n, bn, mas, vel, e)
    integer, intent(in) :: n, bn
    real, intent(in)    :: mas(n), vel(n)
    real, intent(out)   :: e
    integer             :: i

    e = 0.
    do i = bn, n - bn
      e = e + 0.5 * mas(i) * vel(i) ** 2
    end do
  end subroutine get_kinetic_energy

  subroutine get_slength(n, bn, mas, den, slen, sk)
    integer, intent(in) :: n, bn
    real, intent(in)    :: sk, mas(n), den(n)
    real, intent(out)   :: slen(n)
    integer             :: i

    do i = bn, n - bn
      slen(i) = sk * (mas(i) / den(i))
    end do
    call set_periodic(n, bn, slen)
  end subroutine get_slength

  subroutine derivs(n, bn, pos, mas, den, slen, pres, acc, sos, sk)
    integer, intent(in) :: n, bn
    real, intent(in)    :: sos, sk
    real, intent(out)   :: pos(n), mas(n), den(n), slen(n), pres(n), acc(n)
    integer             :: i

    do i = 1, 3
      call get_density(n, bn, pos, mas, den, slen)
      call get_slength(n, bn, mas, den, slen, sk)
    end do
    call equation_of_state(n, bn, den, pres, sos)
    call get_accel(n, bn, pos, mas, den, slen, pres, acc)
  end subroutine derivs

end module internal
