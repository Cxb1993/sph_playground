module arrayresize
  implicit none

  public :: resize

  private

  interface resize
    module procedure i4resize, i8resize, cresize, rresize, resize3r
  end interface

contains
  subroutine i4resize(array, oldsize, newsize)
  ! pure subroutine i4resize(array, oldsize, newsize)
    integer, intent(in)                 :: newsize, oldsize
    integer, intent(inout), allocatable :: array(:)
    integer, allocatable                :: tmp(:)
    integer                             :: i, iterlimit

    iterlimit = min(newsize, oldsize)
    allocate(tmp(newsize))
    ! print*, 'current size: ', size(array), ' | coppied elements: ', iterlimit
    do i=1,iterlimit
      tmp(i) = array(i)
    end do
    deallocate(array)
    allocate(array(newsize))
    array(:) = 0
    do i=1,iterlimit
      array(i) = tmp(i)
    end do
    ! deallocate(tmp)
  end subroutine

  pure subroutine i8resize(array, oldsize, newsize)
    integer, intent(in)                    :: newsize, oldsize
    integer(8), intent(inout), allocatable :: array(:)
    integer(8), allocatable                :: tmp(:)
    integer                                :: i, iterlimit

    iterlimit = min(newsize, oldsize)
    allocate(tmp(newsize))
    do i=1,iterlimit
      tmp(i) = array(i)
    end do
    deallocate(array)
    allocate(array(newsize))
    do i=1,iterlimit
      array(i) = tmp(i)
    end do
    ! deallocate(tmp)
  end subroutine

  pure subroutine rresize(array, oldsize, newsize)
    integer, intent(in)              :: newsize, oldsize
    real, intent(inout), allocatable :: array(:)
    real, allocatable                :: tmp(:)
    integer                          :: i, iterlimit

    iterlimit = min(newsize, oldsize)
    allocate(tmp(newsize))
    do i=1,iterlimit
      tmp(i) = array(i)
    end do
    deallocate(array)
    allocate(array(newsize))
    do i=1,iterlimit
      array(i) = tmp(i)
    end do
    ! deallocate(tmp)
  end subroutine

  pure subroutine cresize(array, chsz, oldsize, newsize)
    integer, intent(in)                          :: newsize, oldsize, chsz
    character(len=*), allocatable, intent(inout) :: array(:)
    character(len=chsz), allocatable             :: tmp(:)
    integer                                      :: i, iterlimit

    iterlimit = min(newsize, oldsize)
    allocate(tmp(newsize))
    do i=1,iterlimit
      tmp(i) = array(i)
    end do
    deallocate(array)
    allocate(array(newsize))
    do i=1,iterlimit
      array(i) = tmp(i)
    end do
    ! deallocate(tmp)
  end subroutine

  pure subroutine resize3r(array, oldsize, newsize)
    integer, intent(in)              :: newsize, oldsize
    real, intent(inout), allocatable :: array(:,:)
    real, allocatable                :: tmp(:,:)
    integer                          :: i, iterlimit

    iterlimit = min(newsize, oldsize)
    allocate(tmp(3,newsize))
    do i=1,iterlimit
      tmp(:,i) = array(:,i)
    end do
    deallocate(array)
    allocate(array(3,newsize))
    do i=1,iterlimit
      array(:,i) = tmp(:,i)
    end do
    ! deallocate(tmp)
  end subroutine
end module
