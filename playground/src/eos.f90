module eos
  use omp_lib
  use state,  only: get_tasktype

  public :: eos_adiabatic, eos_isothermal
  private
contains

  subroutine eos_adiabatic(den, u, P, c, gamma)
    real, allocatable, intent(in)    :: den(:), u(:)
    real, allocatable, intent(inout) :: P(:), c(:)
    real, intent(in)    :: gamma
    integer             :: i

    do i = 1, size(den)
      ! if (tasktype == 4) then
      !   P(i) = P(i) * 1
      ! end if
      P(i) = (gamma - 1) * den(i) * u(i)
      c(i) = sqrt(gamma * P(i) / den(i))
    end do
  end subroutine

  subroutine eos_isothermal(den, c, P)
    real, allocatable, intent(in)    :: den(:)
    real,              intent(in)    :: c
    real, allocatable, intent(inout) :: P(:)
    integer             :: i, n

    n = size(P)

    do i = 1, n
      P(i) = den(i) * c * c
    end do
  end subroutine
end module
