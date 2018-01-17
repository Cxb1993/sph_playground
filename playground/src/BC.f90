module BC
  use list
  use const
  use timing,       only: addTime
  use state,        only: getdim
  use arrayresize,  only: resize

  implicit none

  public :: &
    setBorders, setRealPartNumber, getRealPartNumber, reflecPeriodicParticles,&
    createPeriodicBorder, clearPeriodicParticles, createFixedBorders,&
    findInsideBorderParticles, getCrossRef, getArtPartNumber, updateFixedToSymmetric

  private
  save
    type(intlist) :: &
      ibx1, ibx2, iby1, iby2, ibz1, ibz2, &
      crossref, fixedreal, fixedbrdr
    integer(8) :: &
      start=0, finish=0

    integer :: realpartnumb, fixedpartnumb, periodpartnumb, initdone=0
    real    :: xmin, xmax, ymin, ymax, zmin, zmax, bordersize
    real    :: boxmax(3), boxmin(3)

contains
  subroutine setRealPartNumber(r)
    integer, intent(in) :: r
    realpartnumb = r
  end subroutine
  pure subroutine getRealPartNumber(r)
    integer, intent(out) :: r
    r = realpartnumb
  end subroutine

  pure subroutine getArtPartNumber(p, f)
    integer, intent(out) :: p, f
    p = periodpartnumb
    f = fixedpartnumb
  end subroutine

  pure function getCrossRef(ini) result(oti)
    integer, intent(in) :: ini
    integer oti
    if (ini <= realpartnumb + fixedpartnumb) then
      oti = ini
    else
      oti = crossref%xe(ini-(realpartnumb + fixedpartnumb))
    end if
  end function getCrossRef

  subroutine setBorders(x1,x2,y1,y2,z1,z2,db)
    real, intent(in) :: &
      x1,x2,y1,y2,z1,z2,db

    bordersize = db
    xmin = x1
    xmax = x2
    ymin = y1
    ymax = y2
    zmin = z1
    zmax = z2
    boxmax(:) = [xmax, ymax, zmax] - bordersize
    boxmin(:) = [xmin, ymin, zmin] + bordersize
  end subroutine

  subroutine clearPeriodicParticles(store)
    real, allocatable, intent(inout) :: store(:,:)

    periodpartnumb = 0
    call crossref%clearfast()
    call ibx1%clearfast()
    call ibx2%clearfast()
    call iby1%clearfast()
    call iby2%clearfast()
    call ibz1%clearfast()
    call ibz2%clearfast()
  end subroutine clearPeriodicParticles

  subroutine reflecPeriodicParticles(store, refdir)
    real, allocatable, intent(inout) :: store(:,:)
    integer, intent(in) :: refdir
    integer :: i, dim
    real :: edgemax(3), edgemin(3), dxmin(3), dxmax(3)

    call system_clock(start)

    edgemax(1) = xmax
    edgemax(2) = ymax
    edgemax(3) = zmax

    edgemin(1) = xmin
    edgemin(2) = ymin
    edgemin(3) = zmin

    call getdim(dim)

    select case(refdir)
    case (ebc_all)
      do i=1,realpartnumb
        dxmin(:) = edgemin(:) - store(es_rx:es_rz, i)
        dxmax(:) = store(es_rx:es_rz, i) - edgemax(:)
        if (dxmax(1) > 0) then
          store(es_rx, i) = xmin + dxmax(1)
        end if
        if (dxmin(1) > 0) then
          store(es_rx, i) = xmax - dxmin(1)
        end if
        if ( dim > 1) then
          if (dxmax(2) > 0) then
            store(es_ry, i) = ymin + dxmax(2)
          end if
          if (dxmin(2) > 0) then
            store(es_ry, i) = ymax - dxmin(2)
          end if
        end if
        if (dim == 3) then
          if (dxmax(3) > 0) then
            store(es_rz, i) = zmin + dxmax(3)
          end if
          if (dxmin(3) > 0) then
            store(es_rz, i) = zmax - dxmin(3)
          end if
        end if
      end do
    case (ebc_x)
      do i=1,realpartnumb
        dxmin(:) = edgemin(:) - store(es_rx:es_rz, i)
        dxmax(:) = store(es_rx:es_rz, i) - edgemax(:)
        if (dxmax(1) > 0) then
          store(es_rx, i) = xmin + dxmax(1)
        end if
        if (dxmin(1) > 0) then
          store(es_rx, i) = xmax - dxmin(1)
        end if
      end do
    case (ebc_y)
      if (dim > 1) then
        do i=1,realpartnumb
          dxmin(:) = edgemin(:) - store(es_rx:es_rz, i)
          dxmax(:) = store(es_rx:es_rz, i) - edgemax(:)
          if (dxmax(2) > 0) then
            store(es_ry, i) = ymin + dxmax(2)
          end if
          if (dxmin(2) > 0) then
            store(es_ry, i) = ymax - dxmin(2)
          end if
        end do
      end if
    case (ebc_z)
      if (dim == 3) then
        do i=1,realpartnumb
          dxmin(:) = edgemin(:) - store(es_rx:es_rz, i)
          dxmax(:) = store(es_rx:es_rz, i) - edgemax(:)
          if (dxmax(3) > 0) then
            store(es_rz, i) = zmin + dxmax(3)
          end if
          if (dxmin(3) > 0) then
            store(es_rz, i) = zmax - dxmin(3)
          end if
        end do
      end if
    case default
      print*, '<!!!> desired side reflection is not defined'
      stop
    end select

    call system_clock(finish)
    call addTime(' bc', finish - start)
  end subroutine reflecPeriodicParticles

  subroutine findInsideBorderParticles(store)
    real, allocatable, intent(in) :: store(:,:)
    integer :: i, dim, storesize
    real :: ra(3)

    call system_clock(start)
    call getdim(dim)

    ! TODO
    ! check if it is possible to do with KDT
    storesize = size(store,2)
    if (dim == 1) then
      do i=1,storesize
        if ((store(es_type,i) == ept_real).or.(store(es_type,i) == ept_fixed)) then
          ra(:) = store(es_rx:es_rz,i)
          if (ra(1) < boxmin(1)) then
            call ibx1%append(i)
          else if (ra(1) > boxmax(1)) then
            call ibx2%append(i)
          end if
        end if
      end do
    else if (dim == 2) then
      do i=1,storesize
        if ((store(es_type,i) == ept_real).or.(store(es_type,i) == ept_fixed)) then
          ra(:) = store(es_rx:es_rz,i)
          if (ra(1) < boxmin(1)) then
            call ibx1%append(i)
          else if (ra(1) > boxmax(1)) then
            call ibx2%append(i)
          end if
          if (ra(2) < boxmin(2)) then
            call iby1%append(i)
          else if (ra(2) > boxmax(2)) then
            call iby2%append(i)
          end if
        end if
      end do
    else if (dim == 3) then
      do i=1,storesize
        if ((store(es_type,i) == ept_real).or.(store(es_type,i) == ept_fixed)) then
          ra(:) = store(es_rx:es_rz,i)
          if (ra(1) < boxmin(1)) then
            call ibx1%append(i)
          else if (ra(1) > boxmax(1)) then
            call ibx2%append(i)
          end if
          if (ra(2) < boxmin(2)) then
            call iby1%append(i)
          else if (ra(2) > boxmax(2)) then
            call iby2%append(i)
          end if
          if (ra(3) < boxmin(3)) then
            call ibz1%append(i)
          else if (ra(3) > boxmax(3)) then
            call ibz2%append(i)
          end if
        end if
      end do
    end if
    call system_clock(finish)
    call addTime(' bc', finish - start)
  end subroutine

  subroutine createPeriodicBorder(store, targetSide)
    real, allocatable, intent(inout) :: store(:,:)
    integer, intent(in) :: targetSide

    integer :: i, pi, k, dim, storesize, prevpn, selfcrossref

    call system_clock(start)
    call getdim(dim)
    storesize = size(store,2)
    prevpn = realpartnumb + fixedpartnumb + periodpartnumb
    selfcrossref = realpartnumb + fixedpartnumb
    k = prevpn
    if ((targetSide == ebc_all).or.(targetSide == ebc_x)) then
      do i=1,ibx1%llen()
        k = k + 1
        pi = ibx1%xe(i)
        if (storesize < k) then
          call resize(store, storesize, 2*storesize)
          storesize = 2*storesize
          print*, "# <?> store expanded to", storesize, "in BC.Periodic.X1"
        end if
        store(es_rx:es_rz,k) = store(es_rx:es_rz,pi)
        store(es_rx,k) = xmax + (store(es_rx,k) - xmin)
        store(es_type,k) = ept_periodic
        call crossref%append(pi)
        if (dim > 1) then
          if (store(es_ry,k) < boxmin(2)) then
            call iby1%append(k)
          else if (store(es_ry,k) > boxmax(2)) then
            call iby2%append(k)
          end if
          if (dim == 3) then
            if (store(es_rz,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (store(es_rz,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
        end if
      end do
      do i=1,ibx2%llen()
        k = k + 1
        pi = ibx2%xe(i)
        if (storesize < k) then
          call resize(store, storesize, 2*storesize)
          storesize = 2*storesize
          print*, "# <?> store expanded to", storesize, "in BC.Periodic.X2"
        end if
        store(es_rx:es_rz,k) = store(es_rx:es_rz,pi)
        store(es_rx,k) = xmin + (store(es_rx,k) - xmax)
        store(es_type,k) = ept_periodic
        call crossref%append(pi)
        if (dim > 1) then
          if (store(es_ry,k) < boxmin(2)) then
            call iby1%append(k)
          else if (store(es_ry,k) > boxmax(2)) then
            call iby2%append(k)
          end if
          if (dim == 3) then
            if (store(es_rz,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (store(es_rz,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
        end if
      end do
    end if
    if (dim > 1) then
      if ((targetSide == ebc_all).or.(targetSide == ebc_y)) then
        do i=1,iby1%llen()
          k = k + 1
          pi = iby1%xe(i)
          if (storesize < k) then
            call resize(store, storesize, 2*storesize)
            storesize = 2*storesize
            print*, "# <?> store expanded to", storesize, "in BC.Periodic.Y1"
          end if
          store(es_rx:es_rz,k) = store(es_rx:es_rz,pi)
          store(es_ry,k) = ymax + (store(es_ry,k) - ymin)
          store(es_type,k) = ept_periodic
          do while (pi > selfcrossref)
            pi = getCrossRef(pi)
          end do
          call crossref%append(pi)
          if (dim == 3) then
            if (store(es_rz,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (store(es_rz,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
        end do
        do i=1,iby2%llen()
          k = k + 1
          pi = iby2%xe(i)
          if (storesize < k) then
            call resize(store, storesize, 2*storesize)
            storesize = 2*storesize
            print*, "# <?> store expanded to", storesize, "in BC.Periodic.Y2"
          end if
          store(es_rx:es_rz,k) = store(es_rx:es_rz,pi)
          store(es_ry,k) = ymin + (store(es_ry,k) - ymax)
          store(es_type,k) = ept_periodic
          do while (pi > selfcrossref)
            pi = getCrossRef(pi)
          end do
          call crossref%append(pi)
          if (dim == 3) then
            if (store(es_rz,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (store(es_rz,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
        end do
      end if
      if (dim == 3) then
        if ((targetSide == ebc_all).or.(targetSide == ebc_z)) then
          do i=1,ibz1%llen()
            k = k + 1
            pi = ibz1%xe(i)
            if (storesize < k) then
              call resize(store, storesize, 2*storesize)
              storesize = 2*storesize
              print*, "# <?> store expanded to", storesize, "in BC.Periodic.Z1"
            end if
            store(es_rx:es_rz,k) = store(es_rx:es_rz,pi)
            store(es_rz,k) = zmax + (store(es_rz,k) - zmin)
            store(es_type,k) = ept_periodic
            do while (pi > selfcrossref)
              pi = getCrossRef(pi)
            end do
            call crossref%append(pi)
          end do
          do i=1,ibz2%llen()
            k = k + 1
            pi = ibz2%xe(i)
            if (storesize < k) then
              call resize(store, storesize, 2*storesize)
              storesize = 2*storesize
              print*, "# <?> store expanded to", storesize, "in BC.Periodic.Z2"
            end if
            store(es_rx:es_rz,k) = store(es_rx:es_rz,pi)
            store(es_rz,k) = zmin + (store(es_rz,k) - zmax)
            store(es_type,k) = ept_periodic
            do while (pi > selfcrossref)
              pi = getCrossRef(pi)
            end do
            call crossref%append(pi)
          end do
        end if
      end if
    end if
    periodpartnumb = periodpartnumb + (k - prevpn)
    call system_clock(finish)
    call addTime(' bc', finish - start)
  end subroutine createPeriodicBorder

  subroutine createFixedBorders(store, targetSide)
    real, allocatable, intent(inout) :: store(:,:)
    integer, intent(in) :: targetSide

    integer :: &
      i, pi, k, dim, storesize, refidx, prevrfpn

    call findInsideBorderParticles(store)

    call system_clock(start)
    call getdim(dim)
    storesize = size(store,2)
    prevrfpn = realpartnumb + fixedpartnumb
    k = prevrfpn
    if ((targetSide == ebc_all).or.(targetSide == ebc_x)) then
      do i=1,ibx1%llen()
        k = k + 1
        pi = ibx1%xe(i)
        call fixedreal%append(pi)
        call fixedbrdr%append(k)
        if (storesize < k) then
          call resize(store,storesize, 2*storesize)
          storesize = 2*storesize
          print*, "# <?> store expanded to", storesize, "in BC.Fixed.X1"
        end if
        store(es_rx:es_rz,k) = store(es_rx:es_rz,pi)
        store(es_rx,k) = xmin + (xmin - store(es_rx,k))
        store(es_type,k) = ept_fixed
        if (dim > 1) then
          if (store(es_ry,k) < boxmin(2)) then
            call iby1%append(k)
          else if (store(es_ry,k) > boxmax(2)) then
            call iby2%append(k)
          end if
          if (dim == 3) then
            if (store(es_rz,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (store(es_rz,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
        end if
      end do
      do i=1,ibx2%llen()
        k = k + 1
        pi = ibx2%xe(i)
        call fixedreal%append(pi)
        call fixedbrdr%append(k)
        if (storesize < k) then
          call resize(store,storesize, 2*storesize)
          storesize = 2*storesize
          print*, "# <?> store expanded to", storesize, "in BC.Fixed.X2"
        end if
        store(es_rx:es_rz,k) = store(es_rx:es_rz,pi)
        store(es_rx,k) = xmax + (xmax - store(es_rx,k))
        store(es_type,k) = ept_fixed
        if (dim > 1) then
          if (store(es_ry,k) < boxmin(2)) then
            call iby1%append(k)
          else if (store(es_ry,k) > boxmax(2)) then
            call iby2%append(k)
          end if
          if (dim == 3) then
            if (store(es_rz,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (store(es_rz,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
        end if
      end do
    end if
    if (dim > 1) then
      if ((targetSide == ebc_all).or.(targetSide == ebc_y)) then
        do i=1,iby1%llen()
          k = k + 1
          pi = iby1%xe(i)
          call fixedreal%append(pi)
          call fixedbrdr%append(k)
          if (storesize < k) then
            call resize(store,storesize, 2*storesize)
            storesize = 2*storesize
            print*, "# <?> store expanded to", storesize, "in BC.Fixed.Y1"
          end if
          store(es_rx:es_rz,k) = store(es_rx:es_rz,pi)
          store(es_ry,k) = ymin + (ymin - store(es_ry,k))
          store(es_type,k) = ept_fixed
          if (dim == 3) then
            if (store(es_rz,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (store(es_rz,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
        end do
        do i=1,iby2%llen()
          k = k + 1
          pi = iby2%xe(i)
          call fixedreal%append(pi)
          call fixedbrdr%append(k)
          if (storesize < k) then
            call resize(store,storesize, 2*storesize)
            storesize = 2*storesize
            print*, "# <?> store expanded to", storesize, "in BC.Fixed.Y2"
          end if
          store(es_rx:es_rz,k) = store(es_rx:es_rz,pi)
          store(es_ry,k) = ymax + (ymax - store(es_ry,k))
          store(es_type,k) = ept_fixed
          if (dim == 3) then
            if (store(es_rz,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (store(es_rz,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
        end do
      end if
      if (dim == 3) then
        if ((targetSide == ebc_all).or.(targetSide == ebc_z)) then
          do i=1,ibz1%llen()
            k = k + 1
            pi = ibz1%xe(i)
            call fixedreal%append(pi)
            call fixedbrdr%append(k)
            if (storesize < k) then
              call resize(store,storesize, 2*storesize)
              storesize = 2*storesize
              print*, "# <?> store expanded to", storesize, "in BC.Fixed.Z1"
            end if
            store(es_rx:es_rz,k) = store(es_rx:es_rz,pi)
            store(es_rz,k) = zmin + (zmin - store(es_rz,k))
            store(es_type,k) = ept_fixed
          end do
          do i=1,ibz2%llen()
            k = k + 1
            pi = ibz2%xe(i)
            call fixedreal%append(pi)
            call fixedbrdr%append(k)
            if (storesize < k) then
              call resize(store,storesize, 2*storesize)
              storesize = 2*storesize
              print*, "# <?> store expanded to", storesize, "in BC.Fixed.Z2"
            end if
            store(es_rx:es_rz,k) = store(es_rx:es_rz,pi)
            store(es_rz,k) = zmax + (zmax - store(es_rz,k))
            store(es_type,k) = ept_fixed
          end do
        end if
      end if
    end if
    fixedpartnumb = fixedpartnumb + (k - prevrfpn)
    call system_clock(finish)
    call addTime(' bc', finish - start)
  end subroutine createFixedBorders

  subroutine updateFixedToSymmetric(store, indexes)
    real, allocatable, intent(inout) :: store(:,:)
    integer, dimension(:), intent(in) :: indexes

    integer :: &
      i, j, ri, bi

    call system_clock(start)
    do i = 1,fixedreal%llen()
      ri = fixedreal%xe(i)
      bi = fixedbrdr%xe(i)
      do j = 1,size(indexes,1)
        store(indexes(j),bi) = -store(indexes(j),ri)
      end do
    end do
    call system_clock(finish)
    call addTime(' bc', finish - start)
  end subroutine updateFixedToSymmetric
end module BC
