module utils
  implicit none

  public :: iresize, cresize, rresize, resize3r

  private

contains
  subroutine iresize(array, oldsize, newsize)
    integer, intent(in)                 :: newsize, oldsize
    integer, intent(inout), allocatable :: array(:)
    integer, allocatable                :: tmp(:)
    integer                             :: i

    allocate(tmp(newsize))
    do i=1,oldsize
      tmp(i) = array(i)
    end do
    deallocate(array)
    allocate(array(newsize))
    do i=1,oldsize
      array(i) = tmp(i)
    end do
  end subroutine

  subroutine rresize(array, oldsize, newsize)
    integer, intent(in)              :: newsize, oldsize
    real, intent(inout), allocatable :: array(:)
    real, allocatable                :: tmp(:)
    integer                          :: i

    allocate(tmp(newsize))
    do i=1,oldsize
      tmp(i) = array(i)
    end do
    deallocate(array)
    allocate(array(newsize))
    do i=1,oldsize
      array(i) = tmp(i)
    end do
  end subroutine

  subroutine cresize(array, chsz, oldsize, newsize)
    integer, intent(in)                          :: newsize, oldsize, chsz
    character(len=*), allocatable, intent(inout) :: array(:)
    character(len=chsz), allocatable               :: tmp(:)
    integer                                      :: i

    allocate(tmp(newsize))
    do i=1,oldsize
      tmp(i) = array(i)
    end do
    deallocate(array)
    allocate(array(newsize))
    do i=1,oldsize
      array(i) = tmp(i)
    end do
  end subroutine

  subroutine resize3r(array, oldsize, newsize)
    integer, intent(in)              :: newsize, oldsize
    real, intent(inout), allocatable :: array(:,:)
    real, allocatable                :: tmp(:,:)
    integer                          :: i

    allocate(tmp(3,newsize))
    do i=1,oldsize
      tmp(:,i) = array(:,i)
    end do
    deallocate(array)
    allocate(array(3,newsize))
    do i=1,oldsize
      array(:,i) = tmp(:,i)
    end do
  end subroutine
end module utils
