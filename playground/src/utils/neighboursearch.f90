module NeighbourSearch
  use const
  use Timing,       only: addTime
  use kernel,       only: get_krad
  use state,        only: getdim,&
                          getddwtype
  use ArrayResize,  only: resize
  use list,         only: intlist
  use omp_lib
  use kdtree2_module

  implicit none

public  findneighboursKDT, &
        getneighbours, getNeibNumbers, destroy, &
        setStepsize, isInitialized, getNeibListL1, getNeibListL2,&
        xgetneighindex, xgetneighnumber!,&
        ! findneighboursN2once, findneighboursN2, findneighboursN2plusStatic,

private
save
  type(intlist), allocatable :: neighbours(:)
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
    if (.not.allocated(ol)) then
      allocate(ol(size(alllistlv1)))
    else if (size(alllistlv1) /= size(ol)) then
      call resize(ol, 1, size(alllistlv1))
    end if
    ol(:) = alllistlv1(:)
  end subroutine getNeibListL1

  pure subroutine getNeibListL2(ol)
    integer, allocatable, intent(inout) :: ol(:)
    if (.not.allocated(ol)) then
      allocate(ol(size(alllistlv2)))
    else if (size(alllistlv2) /= size(ol)) then
      call resize(ol, 1, size(alllistlv2))
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

  subroutine findneighboursKDT(store)!,level
    use kernel, only: getkernelnn => getneibnumber
    use state,  only: getddwtype, getPartNumber
    use BC,     only: getPeriodPartNumber
    real, allocatable, intent(in)     :: store(:,:)
    ! real, intent(in) :: level
    type(kdtree2), pointer            :: kdtree
    type(kdtree2_result), allocatable :: kdtree_res(:)
    integer :: &
      maxresultnum, nfound, nlsz, i, j, ktp, &
      realpn, periodpn, fixedpn,&
      al1, al2, sn, ln
    real    :: kr

    call system_clock(start)
    call get_krad(kr)
    call getkernelnn(maxresultnum)
    call getddwtype(ktp)
    call getPartNumber(realpn,fixedpn)
    call getPeriodPartNumber(periodpn)

    sn = realpn + fixedpn + periodpn
    kdtree => kdtree2_create(store(es_rx:es_rz,1:sn))
    allocate(kdtree_res(maxresultnum))
    if (.not.allocated(neighbours)) then
      allocate(neighbours(realpn))
    end if

    if (allocated(alllistlv1)) then
      deallocate(alllistlv1)
    end if
    allocate(alllistlv1(realpn))
    ln = realpn
    ! if (allocated(alllistlv2)) then
    !   deallocate(alllistlv2)
    ! end if
    ! if (level == enl_l1) then
    !   allocate(alllistlv1(realpn))
    !   ln = realpn
    ! else if (level == enl_l2) then
    !   allocate(alllistlv1(realpn+fixedpn))
    !   allocate(alllistlv2(realpn+fixedpn))
    !   ln = realpn + fixedpn
    ! end if

    alllistlv1(:) = 0
    !$omp parallel do default(none)&
    !$omp private(i, j, nlsz, nfound, kdtree_res)&
    !$omp shared(store, kr, neighbours, stepsize, ktp, realpn)&
    !$omp shared(alllistlv1, alllistlv2, maxresultnum, kdtree, ln)
    do i=1,ln,stepsize
      alllistlv1(i) = 1
      call kdtree2_r_nearest_around_point(kdtree, i, 0, (kr * store(es_h,i))**2, nfound, maxresultnum, kdtree_res)
      if (maxresultnum < nfound) then
        print*, "# <!> need at least ", nfound, " particles."
        print*, "# <!> current max is ", maxresultnum
        print*, "# <!> pos = [", store(es_rx:es_rz,i), "]; hfac = ", store(es_h,i)
        error stop "# <!> KDTree neibsearch result was truncated."
      end if
      ! if (.not.(allocated(neighbours(i)%list))) then
      !   allocate(neighbours(i)%list(nfound-1))
      ! end if
      ! nlsz = size(neighbours(i)%list)
      ! if (nlsz /= (nfound-1)) then
      !   call resize(neighbours(i)%list, 1, nfound-1)
      ! end if
      call kdtsearchFiller(i, nfound, kdtree_res)
    end do
    !$omp end parallel do

    al1 = 1
    do i = 1,ln
      if (alllistlv1(i) == 1) then
        alllistlv1(al1) = i
        al1 = al1 + 1
      end if
    end do
    al1 = al1 - 1
    call resize(alllistlv1, al1, al1)

    ! al1 = 1
    ! al2 = 1
    ! do i = 1,ln
    !   if (alllistlv1(i) == 1) then
    !     if (store(es_type) == ept_real) then
    !       alllistlv1(al1) = i
    !       alllistlv2(al2) = i
    !       al1 = al1 + 1
    !       al2 = al2 + 1
    !     else
    !       alllistlv2(al2) = i
    !       al2 = al2 + 1
    !     end if
    !   end if
    ! end do
    ! al1 = al1 - 1
    ! al2 = al2 - 1
    ! call resize(alllistlv1, al1, al1)
    ! call resize(alllistlv2, al2, al2)

    initialized = 1.
    call kdtree2_destroy(kdtree)
    deallocate(kdtree_res)

    call system_clock(finish)
    call addTime(' neibs', finish - start)
  end subroutine findneighboursKDT

  subroutine kdtsearchFiller(i, nfound, kdtree_res)
    type(kdtree2_result), allocatable, intent(inout) :: kdtree_res(:)
    integer, intent(in) :: i, nfound
    integer :: j, k

    call neighbours(i)%clearfast()
    do j = 1,nfound
      k = kdtree_res(j)%idx
      if (k /= i) then
        call neighbours(i)%append(k)
      end if
    end do
  end subroutine kdtsearchFiller

  subroutine getneighbours(idx, list, dt)
    integer(8), intent(inout)           :: dt
    integer, allocatable, intent(inout) :: list(:)
    integer, intent(in)                 :: idx
    integer                             :: snl
    call system_clock(start)

    if (allocated(neighbours)) then
      snl = neighbours(idx)%llen()
      if (allocated(list)) then
        if (size(list) /= snl) then
          call resize(list, size(list), snl)
        end if
      else
        allocate(list(snl))
      end if
      list(:) = neighbours(idx)%toarr()
    else
      error stop "<!> there were no niebs search before"
    end if
    call system_clock(finish)
    dt = finish - start
    call addTime(' neibs', dt)
  end subroutine getneighbours

  pure function xgetneighnumber(ida) result(num)
    integer, intent(in) :: &
      ida
    integer :: &
      num

      num = neighbours(ida)%llen()
  end function xgetneighnumber

  pure function xgetneighindex(ida, idb) result(idx)
    integer, intent(in) :: &
      ida, idb
    integer :: &
      idx

      idx = neighbours(ida)%xe(idb)
  end function xgetneighindex

  subroutine destroy()
    integer :: i, sn

    if(allocated(neighbours)) then
      sn = size(neighbours)
      do i=1,sn
        ! if (allocated(neighbours(i))) then
        call neighbours(i)%clear()
        ! deallocate(neighbours(i))
        ! end if
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
