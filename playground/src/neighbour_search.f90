module neighboursearch
  use kernel, only: get_krad,&
                    get_dim,&
                    get_kerntype
  use utils,  only: resize
  use omp_lib

  implicit none

public findneighbours, getneighbours, setStepsize, isInitialized

private

  type neighbourlisttype
    integer, allocatable :: list(:)
  end type
  type(neighbourlisttype), allocatable :: neighbours(:)

  save
  integer :: stepsize = 1
  integer :: initialized = 0
contains

  subroutine setStepsize(i)
    integer, intent(in) :: i
    stepsize = i
  end subroutine setStepsize

  subroutine isInitialized(o)
    integer, intent(out) :: o
    o = initialized
  end subroutine isInitialized
  ! simple list
  subroutine findneighbours(ptype, pos, h)
    real, allocatable, intent(in)    :: pos(:,:), h(:)
    integer, allocatable, intent(in) :: ptype(:)
    integer, allocatable             :: tmp(:)
    integer                          :: sn, i, j, tsz, tix, dim, kt
    real                             :: r2, r(3), kr

    sn = size(pos, dim=2)
    call get_krad(kr)
    call get_dim(dim)
    call get_kerntype(kt)

    if(allocated(neighbours)) then
      do i=1,sn,stepsize
        deallocate(neighbours(i)%list)
      end do
      deallocate(neighbours)
    end if
    allocate(neighbours(sn))

    !$omp parallel do default(none)&
    !$omp shared(pos, ptype, h, sn, kr, neighbours, dim, stepsize, kt)&
    !$omp private(i, j, tix, r, r2, tsz, tmp)
    do i=1,sn,stepsize
      ! print*, i
      if ((ptype(i) /= 0) .or. (kt == 3)) then
        if (dim == 1) then
          allocate(neighbours(i)%list(10))
        else if (dim == 2) then
          allocate(neighbours(i)%list(50))
        else
          allocate(neighbours(i)%list(100))
        end if
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
        ! print*, tix
        ! read *
      end if
    end do
    !$omp end parallel do
    initialized = 1.
  end subroutine findneighbours

  subroutine getneighbours(idx, pos, h, list)
    real, allocatable, intent(in)       :: pos(:,:), h(:)
    integer, allocatable, intent(out) :: list(:)
    integer, intent(in)               :: idx

    if (allocated(neighbours(idx)%list)) then
      ! print*, 'Try to used old list', idx
      allocate(list(size(neighbours(idx)%list(:))))
      list(:) = neighbours(idx)%list(:)
      ! print*, 'Used old list', idx
    else
      ! print*, 'Try to added new list', idx
      call findneighboursonce(idx, pos, h, list)
      ! print*, 'Added new list', idx
    end if
    ! read*
  end subroutine getneighbours

  subroutine findneighboursonce(idx, pos, h, nlist)
    real, allocatable, intent(in)       :: pos(:,:), h(:)
    integer, allocatable, intent(inout) :: nlist(:)
    integer, intent(in)                 :: idx
    integer                             :: sn, i, j, tsz, tix, dim
    real                                :: r2, r(3), kr

    sn = size(pos, dim=2)
    call get_krad(kr)
    call get_dim(dim)

    if (dim == 1) then
      allocate(nlist(10))
    else if (dim == 2) then
      allocate(nlist(50))
    else
      allocate(nlist(100))
    end if
    tix = 0

    do j = 1,sn
      if ( j /= idx ) then
        r(:) = pos(:,idx) - pos(:,j)
        r2 = dot_product(r(:),r(:))
        if (r2 < (kr * h(idx))**2) then
          tix = tix + 1
          tsz = size(nlist)
          if (tsz < tix) then
            call resize(nlist, tsz, tsz * 2)
          end if
          nlist(tix) = j
        end if
      end if
    end do
    tsz = size(nlist)
    if (tsz /= tix) then
      call resize(nlist, tix, tix)
    end if

    if (.not.allocated(neighbours)) then
      allocate(neighbours(sn))
    end if
    if (.not. allocated(neighbours(idx)%list)) then
      allocate(neighbours(idx)%list(size(nlist)))
    end if
    neighbours(idx)%list(:) = nlist(:)
  end subroutine findneighboursonce
end module neighboursearch
