module utils
  implicit none

  public :: resize

  private

  interface resize
    module procedure i4resize, i8resize, cresize, rresize, resize3r
  end interface

contains
  pure subroutine i4resize(array, oldsize, newsize)
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

  pure subroutine i8resize(array, oldsize, newsize)
    integer, intent(in)                    :: newsize, oldsize
    integer(8), intent(inout), allocatable :: array(:)
    integer(8), allocatable                :: tmp(:)
    integer                                :: i

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

  pure subroutine rresize(array, oldsize, newsize)
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

  pure subroutine cresize(array, chsz, oldsize, newsize)
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

  pure subroutine resize3r(array, oldsize, newsize)
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
