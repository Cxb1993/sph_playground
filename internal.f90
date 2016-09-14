module internal
  use printer
 implicit none

 public :: get_kinetic_energy, derivs

 private

contains

  subroutine set_periodic(n, bn, A)
    integer, intent(in) :: n, bn
    real, intent(out)   :: A(n)
    integer             :: i, nr

    nr = n - 2 * bn
    do i = 1, bn
      A(i) = A(nr + i)
      A(nr + bn + i) = A(bn + i)
    end do
  end subroutine set_periodic

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

  subroutine get_kernel(q, w, dw)
    real, intent(in) :: q
    real             :: w, dw

    if (q >= 2) then
      w  = 0.
      dw = 0.
    else if (q >= 1) then
      w  = 0.25 * ((2. - q) ** 3)
      dw = (- 0.75 * ((2. - q) ** 2)) / q
    else if (q >= 0) then
      w  = 1. - 1.5 * q**2 + 0.75 * q ** 3
      dw = -3. + 2.25 * q
      ! dw = (- 0.75 * ((2. - q) ** 2) + 3.0 * ((1. - q) ** 2)) / q
    else
      print *, 'something went wrong, q =', q
      w  = 0.
      dw = 0.
    end if
     w = 2./3. * w
    dw = 2./3. * dw
  end subroutine get_kernel

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
    real                :: Pi, Pj, wi, wj, dwi, dwj, qi, qj, dx, hi21
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
            dwi = dwi * (pos(i) - pos(j)) / (slen(i) * slen(i))
            dwj = dwj * (pos(i) - pos(j)) / (slen(j) * slen(j))
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
    real, intent(in)    :: sos
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
