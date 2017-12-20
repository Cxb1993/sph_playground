module eos
  use omp_lib
  use state,  only: get_tasktype
  use neighboursearch,  only: getNeibListL1


  public :: eos_adiabatic, eos_isothermal
  private
contains

  subroutine eos_adiabatic(den, u, P, c, gamma)
    real, allocatable, intent(in)     :: den(:), u(:)
    real, allocatable, intent(inout)  :: P(:), c(:)
    real, intent(in)                  :: gamma
    integer, allocatable  :: nlista(:)
    integer               :: la, i

    call getNeibListL1(nlista)
    do la = 1, size(nlista)
      i = nlista(la)
      P(i) = (gamma - 1) * den(i) * u(i)
      c(i) = sqrt(gamma * P(i) / den(i))
    end do
  end subroutine

  subroutine eos_isothermal(den, c, P)
    real, allocatable, intent(in)    :: den(:)
    real,              intent(in)    :: c
    real, allocatable, intent(inout) :: P(:)
    integer, allocatable  :: nlista(:)
    integer               :: la, i

    call getNeibListL1(nlista)
    do la = 1, size(nlista)
      i = nlista(la)
      P(i) = den(i) * c * c
    end do
  end subroutine
end module
