module arrayresize
  implicit none

  public :: resize!, realocnum

  private

  save

  interface resize
    module procedure rsz_i4, i8resize, cresize, rresize, resize3r
  end interface

contains
  pure subroutine rsz_i4(array, oldsize, newsize)
  ! subroutine rsz_i4(array, oldsize, newsize)
    integer, intent(in)                 :: newsize, oldsize
    integer, intent(inout), allocatable :: array(:)
    integer, allocatable                :: tmp(:)
    integer                             :: i, iterlimit

    iterlimit = min(newsize, oldsize)
    allocate(tmp(newsize))
    do i=1,iterlimit
      tmp(i) = array(i)
    end do
    deallocate(array)
    allocate(array(newsize))
    array(iterlimit:newsize) = 0
    do i=1,iterlimit
      array(i) = tmp(i)
    end do
    deallocate(tmp)
    ! realoc = realoc + 1
  end subroutine rsz_i4

  pure subroutine i8resize(array, oldsize, newsize)
  ! subroutine i8resize(array, oldsize, newsize)
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
    array(iterlimit:newsize) = 0
    do i=1,iterlimit
      array(i) = tmp(i)
    end do
    deallocate(tmp)
    ! realoc = realoc + 1
  end subroutine

  pure subroutine rresize(array, oldsize, newsize)
  ! subroutine rresize(array, oldsize, newsize)
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
    array(iterlimit:newsize) = 0.
    do i=1,iterlimit
      array(i) = tmp(i)
    end do
    deallocate(tmp)
    ! realoc = realoc + 1
  end subroutine

  pure subroutine cresize(array, chsz, oldsize, newsize)
  ! subroutine cresize(array, chsz, oldsize, newsize)
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
    array(iterlimit:newsize) = ""
    do i=1,iterlimit
      array(i) = tmp(i)
    end do
    deallocate(tmp)
    ! realoc = realoc + 1
  end subroutine

  pure subroutine resize3r(array, oldsize, newsize)
  ! subroutine resize3r(array, oldsize, newsize)
    integer, intent(in)              :: newsize, oldsize
    real, intent(inout), allocatable :: array(:,:)
    real, allocatable                :: tmp(:,:)
    integer                          :: i, iterlimit, dim1

    dim1 = size(array,dim=1)
    iterlimit = min(newsize, oldsize)
    allocate(tmp(dim1,newsize))
    do i=1,iterlimit
      tmp(:,i) = array(:,i)
    end do
    deallocate(array)
    allocate(array(dim1,newsize))
    array(:,iterlimit:newsize) = 0.
    do i=1,iterlimit
      array(:,i) = tmp(:,i)
    end do
    deallocate(tmp)
    ! realoc = realoc + 1
  end subroutine
end module
