module NeighbourSearch
  use const
  use Timing, only: addTime
  use kernel, only: get_krad
  use state,  only: getdim,&
                    getkerntype
  use ArrayResize,  only: resize
  use omp_lib
  use kdtree2_module

  implicit none

public findneighboursN2, findneighboursN2plusStatic, findneighboursKDT, &
        getneighbours, getNeibNumbers, destroy, &
        setStepsize, isInitialized, getNeibListL1, getNeibListL2, &
        FindPeriodicBorder

private
save
  type neighbourlisttype
    integer, allocatable :: list(:)
  end type
  type(neighbourlisttype), allocatable :: neighbours(:)
  integer, allocatable :: alllistlv1(:), alllistlv2(:)

  integer                :: stepsize = 1
  integer                :: initialized = 0
  integer(8)             :: start = 0, finish = 0
  real, parameter        :: eps = eps0
contains

  subroutine setStepsize(i)
    integer, intent(in) :: i
    stepsize = i
  end subroutine setStepsize

  subroutine isInitialized(o)
    integer, intent(out) :: o
    o = initialized
  end subroutine isInitialized

  pure subroutine getNeibListL1(ol)
    integer, allocatable, intent(inout) :: ol(:)
    if ( .not.allocated(ol) ) then
      allocate(ol(size(alllistlv1)))
    end if
    ol(:) = alllistlv1(:)
  end subroutine getNeibListL1

  pure subroutine getNeibListL2(ol)
    integer, allocatable, intent(inout) :: ol(:)
    if ( .not.allocated(ol) ) then
      allocate(ol(size(alllistlv2)))
    end if
    ol(:) = alllistlv2(:)
  end subroutine getNeibListL2

  pure subroutine getNeibNumbers(ol1, ol2)
    integer, intent(out) :: ol1, ol2

    ol1 = 0
    ol2 = 0
    if (allocated(alllistlv1)) then
      ol1 = size(alllistlv1)
    end if
    if (allocated(alllistlv2)) then
      ol2 = size(alllistlv2)
    end if
  end subroutine

  subroutine FindPeriodicBorder(pos, pspc, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, &
                                bx1, bx2, by1, by2, bz1, bz2, &
                                inbx1, inbx2, inby1, inby2, inbz1, inbz2)
    use list

    real, allocatable, intent(in)     :: pos(:,:)
    real, intent(in)                  :: pspc, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2
    type(intlist), intent(in)         :: bx1, bx2, by1, by2, bz1, bz2
    integer, allocatable, intent(out) :: inbx1(:), inbx2(:), inby1(:), inby2(:), inbz1(:), inbz2(:)

    real    :: approxPos, approxPosVec(3)
    integer :: i, nfound
    type(kdtree2), pointer :: kdtreeStore
    type(kdtree2_result), allocatable :: kdtree_res(:)

    kdtreeStore => kdtree2_create(pos)

    allocate(kdtree_res(10))
    allocate(inbx1(bx1%llen()))
    allocate(inbx2(bx2%llen()))
    allocate(inby1(by1%llen()))
    allocate(inby2(by2%llen()))
    allocate(inbz1(bz1%llen()))
    allocate(inbz2(bz2%llen()))

    do i=1,bx1%llen()
      approxPos = brdx2 - (abs(brdx1 - pos(1,bx1%xe(i))))
      approxPosVec(1) = approxPos
      approxPosVec(2:3) = pos(2:3,bx1%xe(i))
      call kdtree2_r_nearest(kdtreeStore, approxPosVec(:), &
        (pspc/2.)**2, nfound, 10, kdtree_res)
      if (nfound /= 1) then
        error stop "FindPeriodicBorder: Found wrond number of candidates for a periodic particle. FIX IT."
      end if
      inbx1(i) = kdtree_res(1)%idx
      ! print*, '--------------------------'
      ! print*, nfound, kdtree_res(1)%idx, kdtree_res(1)%dis
      ! print*, 'my pos=', pos(:,bz2%xe(i))
      ! print*, 'target pos=', approxPosVec
      ! print*, 'found pos=', pos(:,kdtree_res(1)%idx)
      ! read*
    end do
    do i=1,bx2%llen()
      approxPos = brdx1 + (abs(brdx2 - pos(1,bx2%xe(i))))
      approxPosVec(1) = approxPos
      approxPosVec(2:3) = pos(2:3,bx2%xe(i))
      call kdtree2_r_nearest(kdtreeStore, approxPosVec(:), &
        (pspc/2.)**2, nfound, 10, kdtree_res)
      if (nfound /= 1) then
        error stop "FindPeriodicBorder: Found wrond number of candidates for a periodic particle. FIX IT."
      end if
      inbx2(i) = kdtree_res(1)%idx
    end do

    do i=1,by1%llen()
      approxPos = brdy2 - (abs(brdy1 - pos(2,by1%xe(i))))
      approxPosVec(1) = pos(1,by1%xe(i))
      approxPosVec(2) = approxPos
      approxPosVec(3) = pos(3,by1%xe(i))
      call kdtree2_r_nearest(kdtreeStore, approxPosVec(:), &
        (pspc/2.)**2, nfound, 10, kdtree_res)
      if (nfound /= 1) then
        error stop "FindPeriodicBorder: Found wrond number of candidates for a periodic particle. FIX IT."
      end if
      inby1(i) = kdtree_res(1)%idx
    end do
    do i=1,by2%llen()
      approxPos = brdy1 + (abs(brdy2 - pos(2,by2%xe(i))))
      approxPosVec(1) = pos(1,by2%xe(i))
      approxPosVec(2) = approxPos
      approxPosVec(3) = pos(3,by2%xe(i))
      call kdtree2_r_nearest(kdtreeStore, approxPosVec(:), &
        (pspc/2.)**2, nfound, 10, kdtree_res)
      if (nfound /= 1) then
        error stop "FindPeriodicBorder: Found wrond number of candidates for a periodic particle. FIX IT."
      end if
      inby2(i) = kdtree_res(1)%idx
    end do

    do i=1,bz1%llen()
      approxPos = brdz2 - (abs(brdz1 - pos(3,bz1%xe(i))))
      approxPosVec(1:2) = pos(1:2,bz1%xe(i))
      approxPosVec(3)   = approxPos
      call kdtree2_r_nearest(kdtreeStore, approxPosVec(:), &
        (pspc/2.)**2, nfound, 10, kdtree_res)
      if (nfound /= 1) then
        error stop "FindPeriodicBorder: Found wrond number of candidates for a periodic particle. FIX IT."
      end if
      inbz1(i) = kdtree_res(1)%idx
    end do
    do i=1,bz2%llen()
      approxPos = brdz1 + (abs(brdz2 - pos(3,bz2%xe(i))))
      approxPosVec(1:2) = pos(1:2,bz2%xe(i))
      approxPosVec(3)   = approxPos
      call kdtree2_r_nearest(kdtreeStore, approxPosVec(:), &
        (pspc/2.)**2, nfound, 10, kdtree_res)
      if (nfound /= 1) then
        error stop "FindPeriodicBorder: Found wrond number of candidates for a periodic particle. FIX IT."
      end if
      inbz2(i) = kdtree_res(1)%idx
    end do

    call kdtree2_destroy(kdtreeStore)
  end subroutine FindPeriodicBorder

  subroutine findneighboursKDT(ptype, pos, h)
    use kernel, only: getkernelnn => getneibnumber
    use state,  only: getkerntype

    real, allocatable, intent(in)     :: pos(:,:), h(:)
    integer, allocatable, intent(in)  :: ptype(:)
    type(kdtree2), pointer            :: kdtree
    type(kdtree2_result), allocatable :: kdtree_res(:)
    integer :: maxresultnum, sn, nfound, nlsz, i, j, al1, al2, tx, ktp
    real    :: kr
    call system_clock(start)
    call get_krad(kr)
    call getkernelnn(maxresultnum)
    call getkerntype(ktp)
    ! kdtree2_destroy(kdtree)
    kdtree => kdtree2_create(pos)
    sn = size(pos, dim=2)
    allocate(kdtree_res(maxresultnum))
    if (.not.allocated(neighbours)) then
      allocate(neighbours(sn))
    end if

    if (allocated(alllistlv1)) then
      deallocate(alllistlv1)
      deallocate(alllistlv2)
    end if
    allocate(alllistlv1(sn))
    allocate(alllistlv2(sn))

    alllistlv1(:) = 0
    alllistlv2(:) = 0
    !$omp parallel do default(none)&
    !$omp shared(pos, ptype, h, kr, neighbours, stepsize, ktp)&
    !$omp shared(alllistlv1, alllistlv2, maxresultnum, kdtree)&
    !$omp private(i, j, tx, nlsz, nfound, kdtree_res)
    do i=1,sn,stepsize
      if (ptype(i) /= 0) then
        alllistlv1(i) = 1
        alllistlv2(i) = 1
        tx = 1
        call kdtree2_r_nearest_around_point(kdtree, i, 0, (kr * h(i))**2, nfound, maxresultnum, kdtree_res)
        if (maxresultnum < nfound) then
          print*, "Need at least ", nfound, " particles."
          error stop "KDTree neibsearch result was truncated."
        end if
        if (.not.(allocated(neighbours(i)%list))) then
          allocate(neighbours(i)%list(nfound-1))
        end if
        nlsz = size(neighbours(i)%list)
        if (nlsz /= (nfound-1)) then
          call resize(neighbours(i)%list, 1, nfound-1)
        end if
        do j = 1,nfound
          if (kdtree_res(j)%idx /= i) then
            neighbours(i)%list(tx) = kdtree_res(j)%idx
            alllistlv2(kdtree_res(j)%idx) = 1
            tx = tx + 1
          end if
        end do
      end if
    end do
    !$omp end parallel do

    al1 = 1
    al2 = 1
    do i = 1,sn
      if ( alllistlv1(i) == 1 ) then
        alllistlv1(al1) = i
        al1 = al1 + 1
      end if
      if ( alllistlv2(i) == 1 ) then
        alllistlv2(al2) = i
        al2 = al2 + 1
      end if
    end do
    al1 = al1 - 1
    call resize(alllistlv1, al1, al1)
    al2 = al2 - 1
    call resize(alllistlv2, al2, al2)

    initialized = 1.
    call kdtree2_destroy(kdtree)
    deallocate(kdtree_res)

    call system_clock(finish)
    call addTime(' neibs', finish - start)
  end subroutine

  ! simple list
  subroutine findneighboursN2(ptype, pos, h)
    real, allocatable, intent(in)    :: pos(:,:), h(:)
    integer, allocatable, intent(in) :: ptype(:)
    integer                          :: sn, i, j, tsz, tix, dim, kt, al1, al2
    real                             :: r2, r(3), kr

    call system_clock(start)

    sn = size(pos, dim=2)

    call get_krad(kr)
    call getdim(dim)
    call getkerntype(kt)

    if(allocated(neighbours)) then
      do i=1,sn,stepsize
        if (ptype(i) /= 0) then
          if (allocated(neighbours(i)%list)) then
            deallocate(neighbours(i)%list)
          end if
        end if
      end do
      deallocate(neighbours)
    end if

    allocate(neighbours(sn))

    if (allocated(alllistlv1)) then
      deallocate(alllistlv1)
      deallocate(alllistlv2)
    end if

    allocate(alllistlv1(sn))
    allocate(alllistlv2(sn))

    ! ---------------------
    ! what if ss == 1?
    ! it broke hydroshock
    ! don't remember the reason
    !  to check it
    ! ---------------------
    ! if (stepsize /= 1) then
    alllistlv1(:) = 0
    alllistlv2(:) = 0
    al1 = 1
    al2 = 1
    ! end if

    !$omp parallel do default(none)&
    !$omp shared(pos, ptype, h, sn, kr, neighbours, dim, stepsize, kt)&
    !$omp shared(al1, al2, alllistlv1, alllistlv2)&
    !$omp private(i, j, tix, r, r2, tsz)
    do i=1,sn,stepsize
      if (ptype(i) /= 0) then
        alllistlv2(i) = 1
        if (.not.(allocated(neighbours(i)%list))) then
          if (dim == 1) then
            allocate(neighbours(i)%list(10))
          else if (dim == 2) then
            allocate(neighbours(i)%list(50))
          else
            allocate(neighbours(i)%list(100))
          end if
        end if
        tix = 0
        do j=1,sn
          if (i /= j) then
            r(:) = pos(:,i) - pos(:,j)
            r2 = dot_product(r(:),r(:))
            if (r2 < (kr * h(i))**2 + eps) then
              tix = tix + 1
              tsz = size(neighbours(i)%list)
              if (tsz < tix) then
                call resize(neighbours(i)%list, tsz, tsz * 2)
              end if
              neighbours(i)%list(tix) = j
              alllistlv2(j) = 1
            end if
          end if
        end do
        tsz = size(neighbours(i)%list)
        if (tsz /= tix) then
          call resize(neighbours(i)%list, tix, tix)
        end if
        !$omp critical
        alllistlv1(al1) = i
        al1 = al1 + 1
        !$omp end critical
      end if
    end do
    !$omp end parallel do
    al1 = al1 - 1
    call resize(alllistlv1, al1, al1)
    do i = 1,sn
      if ( alllistlv2(i) == 1 ) then
        alllistlv2(al2) = i
        al2 = al2 + 1
      end if
    end do
    al2 = al2 - 1
    call resize(alllistlv2, al2, al2)

    initialized = 1.

    call system_clock(finish)
    call addTime(' neibs', finish - start)
  end subroutine

  subroutine checkConsistentNeighbours(ptype, pos, h, needDeepCheck)
    real, allocatable, intent(in)       :: pos(:,:), h(:)
    integer, allocatable, intent(in)    :: ptype(:)
    integer, allocatable, intent(inout) :: needDeepCheck(:)
    integer                             :: sn, i, j
    real                                :: r2, r(3), kr

    call get_krad(kr)

    if(.not. allocated(neighbours)) then
      needDeepCheck(:) = 1
    else
      needDeepCheck(:) = 0
      sn = size(pos, dim=2)
      !$omp parallel do default(none)&
      !$omp shared(pos, ptype, h, sn, neighbours, stepsize)&
      !$omp shared(kr, needDeepCheck)&
      !$omp private(i, j, r, r2)
      do i = 1,sn,stepsize
        if (ptype(i) /= 0) then
          if (allocated(neighbours(i)%list)) then
            do j = 1,size(neighbours(i)%list)
              r(:) = pos(:,i) - pos(:,neighbours(i)%list(j))
              r2 = dot_product(r(:),r(:))
              if (r2 > ((kr * h(i))*(kr * h(i)) + eps)) then
                needDeepCheck(i) = 1
                exit
              end if
            end do
          else
            needDeepCheck(i) = 1
          end if
        end if
      end do
      !$omp end parallel do
    end if
    ! print*, needDeepCheck
  end subroutine checkConsistentNeighbours

  subroutine findneighboursN2plusStatic(ptype, pos, h)
    real, allocatable, intent(in)    :: pos(:,:), h(:)
    integer, allocatable, intent(in) :: ptype(:)
    integer, allocatable             :: needDeepCheck(:)
    integer                          :: sn, i, j, tsz, tix, dim, kt, al1, al2
    real                             :: r2, r(3), kr
    call system_clock(start)

    sn = size(pos, dim=2)

    allocate(needDeepCheck(sn))
    call checkConsistentNeighbours(ptype, pos, h, needDeepCheck)

    call get_krad(kr)
    call getdim(dim)
    call getkerntype(kt)

    if (.not. allocated(neighbours)) then
      allocate(neighbours(sn))
    end if

    if (allocated(alllistlv1)) then
      deallocate(alllistlv1)
      deallocate(alllistlv2)
    end if
    allocate(alllistlv1(sn))
    allocate(alllistlv2(sn))

    alllistlv1(:) = 0
    alllistlv2(:) = 0
    al1 = 1
    al2 = 1

    !$omp parallel do default(none)&
    !$omp shared(pos, ptype, h, sn, kr, neighbours, dim, stepsize, kt)&
    !$omp shared(al1, al2, alllistlv1, alllistlv2, needDeepCheck)&
    !$omp private(i, j, tix, r, r2, tsz)
    do i=1,sn,stepsize
      if (ptype(i) /= 0) then
        alllistlv1(i) = 1
        alllistlv2(i) = 1
        if (needDeepCheck(i) == 1) then
          if (allocated(neighbours(i)%list)) then
            deallocate(neighbours(i)%list)
          end if
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
              if (sqrt(r2) < (kr*h(i) + eps)) then
                tix = tix + 1
                tsz = size(neighbours(i)%list)
                if (tsz < tix) then
                  call resize(neighbours(i)%list, tsz, tsz * 2)
                end if
                neighbours(i)%list(tix) = j
                alllistlv2(j) = 1
              end if
            end if
          end do
          tsz = size(neighbours(i)%list)
          if (tsz /= tix) then
            call resize(neighbours(i)%list, tix, tix)
          end if
        end if
      end if
    end do
    !$omp end parallel do

    do i = 1,sn
      if ( alllistlv1(i) == 1 ) then
        alllistlv1(al1) = i
        al1 = al1 + 1
      end if
      if ( alllistlv2(i) == 1 ) then
        alllistlv2(al2) = i
        al2 = al2 + 1
      end if
    end do

    al1 = al1 - 1
    call resize(alllistlv1, al1, al1)
    al2 = al2 - 1
    call resize(alllistlv2, al2, al2)

    initialized = 1.
    call system_clock(finish)
    call addTime(' neibs', finish - start)
  end subroutine

  subroutine getneighbours(idx, pos, h, list, dt)
    real, allocatable, intent(in)       :: pos(:,:), h(:)
    integer(8), intent(inout)           :: dt
    integer, allocatable, intent(inout) :: list(:)
    integer, intent(in)                 :: idx
    integer                             :: sn, snl
    call system_clock(start)

    sn = size(pos, dim=2)
    if (.not.allocated(neighbours)) then
      allocate(neighbours(sn))
    end if
    if (allocated(neighbours(idx)%list)) then
      snl = size(neighbours(idx)%list(:))

      if ( allocated(list) ) then
        if (size(list) /= snl) then
          call resize(list, size(list), snl)
        end if
      else
        allocate(list(snl))
      end if
      list(:) = neighbours(idx)%list(:)
    else
      call findneighboursN2once(idx, pos, h, list)
    end if
    call system_clock(finish)
    dt = finish - start
    call addTime(' neibs', dt)
  end subroutine getneighbours

  subroutine findneighboursN2once(idx, pos, h, nlist)
    real, allocatable, intent(in)       :: pos(:,:), h(:)
    integer, allocatable, intent(inout) :: nlist(:)
    integer, intent(in)                 :: idx
    integer                             :: sn, j, tsz, tix, dim
    real                                :: r2, r(3), kr

    sn = size(pos, dim=2)
    call get_krad(kr)
    call getdim(dim)

    if ( allocated(nlist) ) then
      deallocate(nlist)
    end if
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
        if (r2 < (kr * h(idx))**2 + eps) then
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
  end subroutine findneighboursN2once

  subroutine destroy()
    integer :: i, sn

    if(allocated(neighbours)) then
      sn = size(neighbours)
      do i=1,sn
        if (allocated(neighbours(i)%list)) then
          deallocate(neighbours(i)%list)
        end if
      end do
      deallocate(neighbours)
      if (allocated(alllistlv1)) then
        deallocate(alllistlv1)
      end if
      if (allocated(alllistlv2)) then
        deallocate(alllistlv2)
      end if
    end if
  end subroutine
end module
