module eos
  use omp_lib
  use state,  only: get_tasktype

  public :: eos_adiabatic, eos_isothermal
  private
contains

  subroutine eos_adiabatic(n, den, u, P, c, gamma)
    real, allocatable, intent(in)    :: den(:), u(:)
    real, allocatable, intent(inout) :: P(:), c(:)
    real, intent(in)    :: gamma
    integer, intent(in) :: n
    integer             :: i, t

    call get_tasktype(t)

    !$omp parallel do default(none)&
    !$omp shared(P, den, u, c, n, gamma, t)&
    !$omp private(i)
    do i = 1, n
      ! print *, '--------0'
      P(i) = (gamma - 1) * den(i) * u(i)
      if (t == 4) then
        P(i) = P(i) * 1
      end if
      c(i) = sqrt(gamma * P(i) / den(i))
    end do
    !$omp end parallel do
  end subroutine eos_adiabatic

  subroutine eos_isothermal(n, den, P, c)
    integer, intent(in) :: n
    real, intent(in)    :: den(n), c
    real, intent(out)   :: P(n)
    integer             :: i

    !$omp parallel do default(none)&
    !$omp shared(P, den, c, n)&
    !$omp private(i)
    do i = 1, n
      P(i) = den(i) * c * c
    end do
    !$omp end parallel do
  end subroutine eos_isothermal
end module eos
