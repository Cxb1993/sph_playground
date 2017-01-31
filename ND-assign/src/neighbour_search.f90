module neighboursearch
  use kernel, only: get_krad
  use utils,  only: resize
  use omp_lib

  implicit none

public findneighbours, getneighbours

private

  type neighbourlisttype
    integer, allocatable :: list(:)
  end type
  type(neighbourlisttype), allocatable :: neighbours(:)

contains

  ! simple list
  subroutine findneighbours(pos, h)
    real, allocatable, intent(in) :: pos(:,:), h(:)
    integer                          :: sn, i, j, tsz, tix
    real                             :: r2, r(3), kr

    sn = size(pos, dim=2)
    call get_krad(kr)

    if(allocated(neighbours)) then
      do i=1,sn
        deallocate(neighbours(i)%list)
      end do
      deallocate(neighbours)
    end if
    allocate(neighbours(sn))
    do i=1,sn
      allocate(neighbours(i)%list(10))
    end do

    !$omp parallel do default(none) &
    !$omp shared(pos, h, sn, kr, neighbours) &
    !$omp private(i, j, tix, r, r2, tsz)
    do i=1,sn
      tix = 0
      do j=1,sn
        if (i /= j) then
          r(:) = pos(:,i) - pos(:,j)
          r2 = dot_product(r(:),r(:))
          if (r2 < (kr * h(i))**2) then
            tix = tix + 1
            tsz = size(neighbours(i)%list)
            if (tsz < tix) then
              call resize(neighbours(i)%list, tsz, tsz * 2)
            end if
            neighbours(i)%list(tix) = j
          end if
        end if
      end do
      tsz = size(neighbours(i)%list)
      if (tsz /= tix) then
        call resize(neighbours(i)%list, tix, tix)
      end if
    end do
    !$omp end parallel do
  end subroutine findneighbours

  subroutine getneighbours(i, list)
    integer, allocatable, intent(out)   :: list(:)
    integer, intent(in)                 :: i

    allocate(list(size(neighbours(i)%list(:))))
    list(:) = neighbours(i)%list(:)
  end subroutine getneighbours
end module neighboursearch
