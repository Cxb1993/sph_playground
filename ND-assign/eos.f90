module eos
  use omp_lib

  public :: eos_adiabatic, eos_isothermal
  private
contains

  subroutine eos_adiabatic(n, den, u, P, c, gamma)
    integer, intent(in) :: n
    real, intent(in)    :: den(n), u(n), gamma
    real, intent(out)   :: P(n), c(n)
    integer             :: i

    !$OMP PARALLEL
    !$OMP DO
    do i = 1, n
      P(i) = (gamma - 1) * den(i) * u(i)
      c(i) = sqrt(gamma * P(i) / den(i))
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine eos_adiabatic

  subroutine eos_isothermal(n, den, P, c)
    integer, intent(in) :: n
    real, intent(in)    :: den(n), c
    real, intent(out)   :: P(n)
    integer             :: i

    !$OMP PARALLEL
    !$OMP DO
    do i = 1, n
      P(i) = den(i) * c * c
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine eos_isothermal
end module eos
