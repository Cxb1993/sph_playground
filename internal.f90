module internal
 implicit none

 public :: get_density, equation_of_state, get_accel, get_kinetic_energy
 public :: get_slength, derivs
 private

contains

subroutine get_density(pos, mas, den, slen, n, bn)
  integer, intent(in)            :: n, bn
  real, intent(in)   :: pos(0:n), mas(0:n), slen(0:n)
  real, intent(out)  :: den(0:n)
  real               :: w, dw, q
  integer  :: i, j

  do i = 0, n
    den(i) = 0.0
    do j = 0, n
       q = abs(pos(i)-pos(j))/slen(i)
       if (q < 2.) then
          call get_kernel(q, w, dw)
          den(i) = den(i) + mas(j) * w / slen(i)
       endif
    end do
  end do
  do i = 0, bn - 1
    den(bn - 1 - i) = den(n - bn - i)
    den(n - bn + 1 + i) = den(bn + i)
  end do
end subroutine get_density

subroutine get_kernel(q, w, dw)
  real, intent(in) :: q
  real             :: w, dw

  if (q >= 2) then
    w  = 0.0
    dw = 0.0
  else if (q >= 1) then
    w  = 0.25 * ((2. - q) ** 3)
    dw = (- 0.75 * ((2. - q) ** 2)) / q
  else if (q >= 0) then
    w  = 0.25 * ((2. - q) ** 3) - ((1. - q) ** 3)
    dw = (- 0.75 * ((2. - q) ** 2) + 3.0 * ((1. - q) ** 2)) / q
  else
    print *, 'something went wrong, q =', q
    w  = 0.0
    dw = 0.0
  end if
  w  = 2.0/3. * w
  dw = 2.0/3. * dw
end subroutine get_kernel

subroutine equation_of_state(den, P, c, n, bn)
  integer, intent(in)            :: n, bn
  real, intent(in)   :: den(0:n), c
  real, intent(out)  :: P(0:n)
  integer  :: i

  do i = bn, n - bn
    P(i) = den(i) * c * c
    ! write(*, '(*(F28.24))') float(i), P(i), den(i), c
  end do
  do i = 0, bn - 1
    P(bn - 1 - i) = P(n - bn - i)
    P(n - bn + 1 + i) = P(bn + i)
  end do
end subroutine equation_of_state

subroutine get_accel(pos, mas, den, slen, P, acc, n, bn)
  integer, intent(in)            :: n, bn
  real, intent(in)   :: pos(0:n), mas(0:n), slen(0:n), den(0:n), P(0:n)
  real, intent(out)  :: acc(0:n)
  real               :: P1, P2, w1, w2, dw1, dw2, qi, qj, dx
  integer  :: i, j

  do i = bn, n - bn
    acc(i) = 0.0
    do j = 0, n
      if (i /= j) then
        dx = abs(pos(i)-pos(j))
        qi = dx/slen(i)
        qj = dx/slen(j)
        if (qi < 2. .or. qj < 2.) then
          call get_kernel(qi, w1, dw1)
          call get_kernel(qj, w2, dw2)
          dw1 = dw1 / (slen(i) * slen(i))
          dw2 = dw2 / (slen(j) * slen(j))
          P1 = P(i) * dw1 / (den(i) * den(i))
          P2 = P(j) * dw2 / (den(j) * den(j))
          acc(i) = acc(i) - mas(j) * (P1 + P2)
        end if
      end if
    end do
    ! write(*, '(*(A100))') '-------------------'
    ! write(*, *) i, acc(i)
  end do
  do i = 0, bn - 1
    acc(bn - 1 - i) = acc(n - bn - i)
    acc(n - bn + 1 + i) = acc(bn + i)
  end do
end subroutine get_accel

subroutine get_kinetic_energy(mas, vel, e, n, bn)
  integer, intent(in)            :: n, bn
  real, intent(in)   :: mas(0:n), vel(0:n)
  real, intent(out)  :: e
  integer  :: i

  e = 0.0
  do i = bn, n - bn
    e = e + 1.0/2 * mas(i) * vel(i) * vel(i)
  end do
end subroutine get_kinetic_energy

subroutine get_slength(mas, den, slen, n, bn)
  integer, intent(in)            :: n, bn
  real, intent(in)   :: mas(0:n), den(0:n)
  real, intent(out)  :: slen(0:n)
  integer  :: i

  do i = bn, n - bn
    slen(i) = slen(i) * (mas(i) / den(i))
  end do
  do i = 0, bn - 1
    slen(bn - 1 - i) = slen(n - bn - i)
    slen(n - bn + 1 + i) = slen(bn + i)
  end do
end subroutine get_slength

subroutine derivs(pos, mas, den, slen, pres, acc, n, bn, sos)
  integer          :: n, bn
  real :: sos
  real :: pos(0:n), mas(0:n), den(0:n), slen(0:n), pres(0:n), acc(0:n)
  integer :: i
  ! call get_density(pos, mas, den, slen, n, bn)
  ! call get_slength(mas, den, slen, n, bn)
  do i = 1,3
    call get_density(pos, mas, den, slen, n, bn)
    call get_slength(mas, den, slen, n, bn)
  end do

  call equation_of_state(den, pres, sos, n, bn)
  call get_accel(pos, mas, den, slen, pres, acc, n, bn)
end subroutine derivs

end module internal
