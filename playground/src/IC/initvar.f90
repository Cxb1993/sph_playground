module initvar
  use omp_lib

  use state,        only: get_equations,&
                          getddwtype,&
                          getdim
  use BC

  implicit none

  ! public :: setupShock12

  private
    integer(8) :: start=0, finish=0

contains

  ! subroutine setupShock12(from, to, store)
  !
  ! end subroutine setupShock12

end module
