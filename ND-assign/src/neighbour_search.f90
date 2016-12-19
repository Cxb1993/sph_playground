module neighboursearch
  use kernel
  use utils

  implicit none

public findneighbours, getneighbours

private

  type neighbourlisttype
    integer, allocatable :: list(:)
  end type
  type(neighbourlisttype), allocatable :: neighbour(:)

contains

  subroutine findneighbours(pos, h)
    integer, allocatable, intent(in) :: pos(:,:), h(:)
    integer                          :: sn, i, j, tsz, tix
    real                             :: r2, r(3), kr

    sn = size(pos, dim=2)
    call get_krad(kr)

    if(allocated(neighbour)) then
      do i=1,sn
        deallocate(neighbour(i)%list)
      end do
      deallocate(neighbour)
    end if
    allocate(neighbour(sn))
    do i=1,sn
      allocate(neighbour(i)%list(1))
    end do

    do i=1,sn
      tix = 0
      do j=1,sn
        if (i /= j) then
          r(:) = pos(:,i) - pos(:,j)
          r2 = dot_product(r(:),r(:))
          if (r2 < (kr * h(i))**2) then
            tix = tix + 1
            tsz = size(neighbour(i)%list)
            if (tsz < tix) then
              call resize(neighbour(i)%list, tix, tix *2)
            end if
            neighbour(i)%list(tix) = j
          end if
        end if
      end do
      tsz = size(neighbour(i)%list)
      if (tsz /= tix) then
        call resize(neighbour(i)%list, tix, tixs)
      end if
    end do
  end subroutine findneighbours

  subroutine getneighbours(i, list)
    integer, allocatable, intent(out) :: list(:)
    integer, allocatable, intent(in)  :: i

    list = neighbour(i)%list
  end subroutine getneighbours
end module neighboursearch
