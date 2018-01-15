module eos
  use omp_lib
  use const
  use state,  only: get_tasktype
  use neighboursearch,  only: getNeibListL1


  public :: eos_adiabatic, eos_isothermal
  private
contains

  subroutine eos_adiabatic(store, gamma)
    real, allocatable, intent(inout)  :: store(:,:)
    real, intent(in)                  :: gamma
    integer, allocatable  :: nlista(:)
    integer               :: la, i

    call getNeibListL1(nlista)
    do la = 1, size(nlista)
      i = nlista(la)
      store(es_p,i) = (gamma - 1) * store(es_den,i) * store(es_u,i)
      store(es_c,i) = sqrt(gamma * store(es_p,i) / store(es_den,i))
    end do
  end subroutine

  subroutine eos_isothermal(store)
    real, allocatable, intent(inout) :: store(:,:)
    integer, allocatable  :: nlista(:)
    integer               :: la, i

    call getNeibListL1(nlista)
    do la = 1, size(nlista)
      i = nlista(la)
      store(es_p,i) = store(es_den,i) * store(es_c, 1) * store(es_c, 1)
    end do
  end subroutine
end module
