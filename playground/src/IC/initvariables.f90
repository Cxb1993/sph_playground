module InitVariables
  use omp_lib

  use timing,       only: addTime
  use ArrayResize,  only: resize
  use const
  use kernel,       only: get_krad
  use state,        only: get_tasktype,&
                          getddwtype,&
                          getdim
  use BC

  implicit none

  ! public :: setupIC

  private
  ! integer(8) :: start=0, finish=0
contains
end module
