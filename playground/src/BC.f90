module BC
  use list
  use const
  use timing,       only: addTime
  use state,        only: getdim
  use arrayresize,  only: resize

  implicit none

  public :: destroy, getCrossRef, &
            setBorder, setBorderInside, &
            reflecParticlesPeriodic, createPhantomPeriodic, &
            clearPhantomParticles, createPhantomFixed

  public :: realpartnumb, bordersize, artpartnumb, needcrosref, &
            xmin, xmax, ymin, ymax, zmin, zmax

  private
  save
    integer, allocatable  :: bX1ins(:), bY1ins(:), bZ1ins(:), bX2ins(:), bY2ins(:), bZ2ins(:)
    integer, allocatable  :: bX1exc(:), bY1exc(:), bZ1exc(:), bX2exc(:), bY2exc(:), bZ2exc(:)
    type(intlist)         :: ibx1, ibx2, iby1, iby2, ibz1, ibz2
    type(intlist)         :: ebx1, ebx2, eby1, eby2, ebz1, ebz2
    type(intlist)         :: posref
    integer(8)            :: start=0, finish=0
    ! GLOABL =(
    integer :: realpartnumb=0, artpartnumb=0, initdone=0, needcrosref=0, phantomcreated=1
    real    :: xmin, xmax, ymin, ymax, zmin, zmax, bordersize
    real :: boxmax(3), boxmin(3)

contains
  subroutine clearPhantomParticles(store)
    real, allocatable, intent(inout) :: store(:,:)

    artpartnumb = 0
    needcrosref = 0

    call posref%clearfast()

    call ibx1%clearfast()
    call ibx2%clearfast()
    call iby1%clearfast()
    call iby2%clearfast()
    call ibz1%clearfast()
    call ibz2%clearfast()

    call ebx1%clearfast()
    call ebx2%clearfast()
    call eby1%clearfast()
    call eby2%clearfast()
    call ebz1%clearfast()
    call ebz2%clearfast()

    boxmax(1) = xmax
    boxmax(2) = ymax
    boxmax(3) = zmax

    boxmin(1) = xmin
    boxmin(2) = ymin
    boxmin(3) = zmin

    boxmax(:) = boxmax(:) - bordersize
    boxmin(:) = boxmin(:) + bordersize
    call findInsideBorderParticles(store, boxmax, boxmin)
  end subroutine clearPhantomParticles

  subroutine destroy
    if (initdone == 1) then
      deallocate(bX1exc)
      deallocate(bX2exc)
      deallocate(bY1exc)
      deallocate(bY2exc)
      deallocate(bZ1exc)
      deallocate(bZ2exc)

      deallocate(bX1ins)
      deallocate(bX2ins)
      deallocate(bY1ins)
      deallocate(bY2ins)
      deallocate(bZ1ins)
      deallocate(bZ2ins)
      initdone = 0
    end if
  end subroutine

  pure function getCrossRef(ini) result(oti)
    integer, intent(in) :: ini
    integer oti
    oti = posref%xe(ini)
  end function getCrossRef

  subroutine reflecParticlesPeriodic(pos, refdir)
    real, allocatable, intent(inout) :: pos(:,:)
    integer, intent(in) :: refdir
    integer :: i
    real :: boxmax(3), boxmin(3), dxmin(3), dxmax(3)

    boxmax(1) = xmax
    boxmax(2) = ymax
    boxmax(3) = zmax

    boxmin(1) = xmin
    boxmin(2) = ymin
    boxmin(3) = zmin

    select case(refdir)
    case (ebc_all)
      do i=1,realpartnumb
        dxmin(:) = boxmin(:) - pos(:,i)
        dxmax(:) = pos(:,i) - boxmax(:)
        if (dxmax(1) > 0) then
          pos(1,i) = xmin + dxmax(1)
        end if
        if (dxmin(1) > 0) then
          pos(1,i) = xmax - dxmin(1)
        end if
        if (dxmax(2) > 0) then
          pos(2,i) = ymin + dxmax(2)
        end if
        if (dxmin(2) > 0) then
          pos(2,i) = ymax - dxmin(2)
        end if
        if (dxmax(3) > 0) then
          pos(3,i) = zmin + dxmax(3)
        end if
        if (dxmin(3) > 0) then
          pos(3,i) = zmax - dxmin(3)
        end if
      end do
    case (ebc_x)
      do i=1,realpartnumb
        dxmin(:) = boxmin(:) - pos(:,i)
        dxmax(:) = pos(:,i) - boxmax(:)
        if (dxmax(1) > 0) then
          pos(1,i) = xmin + dxmax(1)
        end if
        if (dxmin(1) > 0) then
          pos(1,i) = xmax - dxmin(1)
        end if
      end do
    case default
      print*, '<!!!> desired side reflection is not defined'
      stop
    end select
  end subroutine reflecParticlesPeriodic

  subroutine findInsideBorderParticles(store, boxmax, boxmin)
    real, allocatable, intent(in) :: store(:,:)
    real, intent(in)    :: boxmax(3), boxmin(3)
    integer :: i, dim
    real :: ra(3)

    call getdim(dim)

    if (dim == 1) then
      do i=1,realpartnumb
        ra(:) = store(es_rx:es_rz,i)
        call posref%append(i)
        if (ra(1) < boxmin(1)) then
          call ibx1%append(i)
        else if (ra(1) > boxmax(1)) then
          call ibx2%append(i)
        end if
      end do
    else if (dim == 2) then
      do i=1,realpartnumb
        ra(:) = store(es_rx:es_rz,i)
        call posref%append(i)
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
      end do
    else if (dim == 3) then
      do i=1,realpartnumb
        ra(:) = store(es_rx:es_rz,i)
        call posref%append(i)
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
      end do
    end if
  end subroutine

  subroutine createPhantomPeriodic(store, targetSide)
    real, allocatable, intent(inout) :: store(:,:)
    integer, intent(in) :: targetSide

    integer :: i, k, dim, storesize, refidx

    call getdim(dim)
    storesize = size(store,dim=2)
    artpartnumb = artpartnumb + 1
    if ((targetSide == ebc_all).or.(targetSide == ebc_x)) then
      do i=1,ibx1%llen()
        k = realpartnumb+artpartnumb
        if (storesize < k) then
          print*, "  <?> pos array was expanded due to phantom periodic particles. L:182"
          call resize(store, k-1, 2*k)
          storesize = 2*k
        end if
        store(:,k) = store(:,ibx1%xe(i))
        store(es_rx,k) = xmax + (store(es_rx,k) - xmin)
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
        artpartnumb = artpartnumb + 1
      end do
      do i=1,ibx2%llen()
        k = realpartnumb+artpartnumb
        if (storesize < k) then
          print*, "  <!> pos array was expanded due to phantom periodic particles. L:208"
          call resize(store, k-1, 2*k)
          storesize = 2*k
        end if
        store(:,k) = store(:,ibx2%xe(i))
        store(es_rx,k) = xmin + (store(es_rx,k) - xmax)
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
        artpartnumb = artpartnumb + 1
      end do
    end if
    if (dim > 1) then
      if ((targetSide == ebc_all).or.(targetSide == ebc_y)) then
        do i=1,iby1%llen()
          k = realpartnumb+artpartnumb
          if (storesize < k) then
            print*, "  <!> pos array was expanded due to phantom periodic particles. L:235"
            call resize(store, k-1, 2*k)
            storesize = 2*k
          end if
          store(:,k) = store(:,iby1%xe(i))
          store(es_ry,k) = ymax + (store(es_ry,k) - ymin)
          if (dim == 3) then
            if (store(es_rz,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (store(es_rz,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
          artpartnumb = artpartnumb + 1
        end do
        do i=1,iby2%llen()
          k = realpartnumb+artpartnumb
          if (storesize < k) then
            print*, "  <!> pos array was expanded due to phantom periodic particles. L:260"
            call resize(store, k-1, 2*k)
            storesize = 2*k
          end if
          store(:,k) = store(:,iby2%xe(i))
          store(es_ry,k) = ymin + (store(es_ry,k) - ymax)
          if (dim == 3) then
            if (store(es_rz,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (store(es_rz,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
          artpartnumb = artpartnumb + 1
        end do
      end if
      if (dim == 3) then
        if ((targetSide == ebc_all).or.(targetSide == ebc_z)) then
          do i=1,ibz1%llen()
            k = realpartnumb+artpartnumb
            if (storesize < k) then
              print*, "  <!> pos array was expanded due to phantom periodic particles. L:286"
              call resize(store, k-1, 2*k)
              storesize = 2*k
            end if
            store(:,k) = store(:,ibz1%xe(i))
            store(es_rz,k) = zmax + (store(es_rz,k) - zmin)
            artpartnumb = artpartnumb + 1
          end do
          do i=1,ibz2%llen()
            k = realpartnumb+artpartnumb
            if (storesize < k) then
              print*, "  <!> pos array was expanded due to phantom periodic particles. L:304"
              call resize(store, k-1, 2*k)
              storesize = 2*k
            end if
            store(:,k) = store(:,ibz2%xe(i))
            store(es_rz,k) = zmin + (store(es_rz,k) - zmax)
            artpartnumb = artpartnumb + 1
          end do
        end if
      end if
    end if
    artpartnumb = artpartnumb - 1
  end subroutine createPhantomPeriodic

  subroutine createPhantomFixed(pos, targetSide)
    real, allocatable, intent(inout) :: pos(:,:)
    integer, intent(in) :: targetSide

    integer :: i, k, dim, spos, refidx

    needcrosref = 1
    call getdim(dim)
    spos = size(pos,dim=2)
    artpartnumb = artpartnumb + 1
    if ((targetSide == ebc_all).or.(targetSide == ebc_x)) then
      do i=1,ibx1%llen()
        k = realpartnumb+artpartnumb
        if (spos < k) then
          print*, "  <!> smal arrays. L:382"
          stop
        end if
        call ebx1%append(k)
        call posref%append(ibx1%xe(i))
        ! call posref%append(k)
        pos(:,k) = pos(:,ibx1%xe(i))
        pos(1,k) = xmin + (xmin - pos(1,k))
        if (dim > 1) then
          if (pos(2,k) < boxmin(2)) then
            call iby1%append(k)
          else if (pos(2,k) > boxmax(2)) then
            call iby2%append(k)
          end if
          if (dim == 3) then
            if (pos(3,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (pos(3,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
        end if
        artpartnumb = artpartnumb + 1
      end do
      do i=1,ibx2%llen()
        k = realpartnumb+artpartnumb
        if (spos < k) then
          print*, "  <!> smal arrays. L:415"
          stop
        end if
        call ebx2%append(k)
        call posref%append(ibx2%xe(i))
        ! call posref%append(k)
        pos(:,k) = pos(:,ibx2%xe(i))
        pos(1,k) = xmax + (xmax - pos(1,k))
        if (dim > 1) then
          if (pos(2,k) < boxmin(2)) then
            call iby1%append(k)
          else if (pos(2,k) > boxmax(2)) then
            call iby2%append(k)
          end if
          if (dim == 3) then
            if (pos(3,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (pos(3,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
        end if
        artpartnumb = artpartnumb + 1
      end do
    end if
    if (dim > 1) then
      if ((targetSide == ebc_all).or.(targetSide == ebc_y)) then
        do i=1,iby1%llen()
          k = realpartnumb+artpartnumb
          if (spos < k) then
            print*, "  <!> smal arrays. L:444"
            stop
          end if
          call eby1%append(k)
          refidx = iby1%xe(i)
          do while (refidx > realpartnumb)
            refidx = getCrossRef(refidx)
          end do
          call posref%append(refidx)
          ! call posref%append(k)
          pos(:,k) = pos(:,iby1%xe(i))
          pos(2,k) = ymin + (ymin - pos(2,k))
          if (dim == 3) then
            if (pos(3,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (pos(3,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
          artpartnumb = artpartnumb + 1
        end do
        do i=1,iby2%llen()
          k = realpartnumb+artpartnumb
          if (spos < k) then
            print*, "  <!> smal arrays. L:463"
            stop
          end if
          call eby2%append(k)
          refidx = iby2%xe(i)
          do while (refidx > realpartnumb)
            refidx = getCrossRef(refidx)
          end do
          call posref%append(refidx)
          ! call posref%append(k)
          pos(:,k) = pos(:,iby2%xe(i))
          pos(2,k) = ymax + (ymax - pos(2,k))
          if (dim == 3) then
            if (pos(3,k) < boxmin(3)) then
              call ibz1%append(k)
            else if (pos(3,k) > boxmax(3)) then
              call ibz2%append(k)
            end if
          end if
          artpartnumb = artpartnumb + 1
        end do
      end if
      if (dim == 3) then
        if ((targetSide == ebc_all).or.(targetSide == ebc_z)) then
          do i=1,ibz1%llen()
            k = realpartnumb+artpartnumb
            if (spos < k) then
              print*, "  <!> smal arrays. L:485"
              stop
            end if
            call ebz1%append(k)
            refidx = ibz1%xe(i)
            do while (refidx > realpartnumb)
              refidx = getCrossRef(refidx)
            end do
            call posref%append(refidx)
            ! call posref%append(k)
            pos(:,k) = pos(:,ibz1%xe(i))
            pos(3,k) = zmin + (zmin - pos(3,k))
            artpartnumb = artpartnumb + 1
          end do
          do i=1,ibz2%llen()
            k = realpartnumb+artpartnumb
            if (spos < k) then
              print*, "  <!> smal arrays. L:497"
              stop
            end if
            call ebz2%append(k)
            refidx = ibz2%xe(i)
            do while (refidx > realpartnumb)
              refidx = getCrossRef(refidx)
            end do
            call posref%append(refidx)
            ! call posref%append(k)
            pos(:,k) = pos(:,ibz2%xe(i))
            pos(3,k) = zmax + (zmax - pos(3,k))
            artpartnumb = artpartnumb + 1
          end do
        end if
      end if
    end if
    artpartnumb = artpartnumb - 1
  end subroutine createPhantomFixed

  ! subroutine initPhantomFixed(targetSide)
  !   integer, intent(in) :: targetSide
  !   integer :: i, k, dim, spos
  !
  !   call getdim(dim)
  !   spos = size(pos,dim=2)
  !   artpartnumb = artpartnumb + 1
  !   if ((targetSide == ebc_all).or.(targetSide == ebc_x)) then
  !     do i=1,ibx1%llen()
  !       k = realpartnumb+artpartnumb
  !       if (spos < k) then
  !         print*, "  <!> smal arrays. L:382"
  !         stop
  !       end if
  !       call ebx1%append(k)
  !       call posref%append(k)
  !       pos(:,k) = pos(:,ibx1%xe(i))
  !       pos(1,k) = xmin + (xmin - pos(1,k))
  !       if (dim > 1) then
  !         if (pos(2,k) < boxmin(2)) then
  !           call iby1%append(k)
  !         else if (pos(2,k) > boxmax(2)) then
  !           call iby2%append(k)
  !         end if
  !         if (dim == 3) then
  !           if (pos(3,k) < boxmin(3)) then
  !             call ibz1%append(k)
  !           else if (pos(3,k) > boxmax(3)) then
  !             call ibz2%append(k)
  !           end if
  !         end if
  !       end if
  !       artpartnumb = artpartnumb + 1
  !     end do
  !     do i=1,ibx2%llen()
  !       k = realpartnumb+artpartnumb
  !       if (spos < k) then
  !         print*, "  <!> smal arrays. L:415"
  !         stop
  !       end if
  !       call ebx2%append(k)
  !       call posref%append(k)
  !       pos(:,k) = pos(:,ibx2%xe(i))
  !       pos(1,k) = xmax + (xmax - pos(1,k))
  !       if (dim > 1) then
  !         if (pos(2,k) < boxmin(2)) then
  !           call iby1%append(k)
  !         else if (pos(2,k) > boxmax(2)) then
  !           call iby2%append(k)
  !         end if
  !         if (dim == 3) then
  !           if (pos(3,k) < boxmin(3)) then
  !             call ibz1%append(k)
  !           else if (pos(3,k) > boxmax(3)) then
  !             call ibz2%append(k)
  !           end if
  !         end if
  !       end if
  !       artpartnumb = artpartnumb + 1
  !     end do
  !   end if
  !   if (dim > 1) then
  !     if ((targetSide == ebc_all).or.(targetSide == ebc_y)) then
  !       do i=1,iby1%llen()
  !         k = realpartnumb+artpartnumb
  !         if (spos < k) then
  !           print*, "  <!> smal arrays. L:444"
  !           stop
  !         end if
  !         call eby1%append(k)
  !         call posref%append(k)
  !         pos(:,k) = pos(:,iby1%xe(i))
  !         pos(2,k) = ymin + (ymin - pos(2,k))
  !         if (dim == 3) then
  !           if (pos(3,k) < boxmin(3)) then
  !             call ibz1%append(k)
  !           else if (pos(3,k) > boxmax(3)) then
  !             call ibz2%append(k)
  !           end if
  !         end if
  !         artpartnumb = artpartnumb + 1
  !       end do
  !       do i=1,iby2%llen()
  !         k = realpartnumb+artpartnumb
  !         if (spos < k) then
  !           print*, "  <!> smal arrays. L:463"
  !           stop
  !         end if
  !         call eby2%append(k)
  !         call posref%append(k)
  !         pos(:,k) = pos(:,iby2%xe(i))
  !         pos(2,k) = ymax + (ymax - pos(2,k))
  !         if (dim == 3) then
  !           if (pos(3,k) < boxmin(3)) then
  !             call ibz1%append(k)
  !           else if (pos(3,k) > boxmax(3)) then
  !             call ibz2%append(k)
  !           end if
  !         end if
  !         artpartnumb = artpartnumb + 1
  !       end do
  !     end if
  !     if (dim == 3) then
  !       if ((targetSide == ebc_all).or.(targetSide == ebc_z)) then
  !         do i=1,ibz1%llen()
  !           k = realpartnumb+artpartnumb
  !           if (spos < k) then
  !             print*, "  <!> smal arrays. L:485"
  !             stop
  !           end if
  !           call ebz1%append(k)
  !           call posref%append(k)
  !           pos(:,k) = pos(:,ibz1%xe(i))
  !           pos(3,k) = zmin + (zmin - pos(3,k))
  !           artpartnumb = artpartnumb + 1
  !         end do
  !         do i=1,ibz2%llen()
  !           k = realpartnumb+artpartnumb
  !           if (spos < k) then
  !             print*, "  <!> smal arrays. L:497"
  !             stop
  !           end if
  !           call ebz2%append(k)
  !           call posref%append(k)
  !           pos(:,k) = pos(:,ibz2%xe(i))
  !           pos(3,k) = zmax + (zmax - pos(3,k))
  !           artpartnumb = artpartnumb + 1
  !         end do
  !       end if
  !     end if
  !   end if
  !   artpartnumb = artpartnumb - 1
  ! end subroutine

  subroutine setBorder(itb, Aexc)
    integer, allocatable, intent(in) :: Aexc(:)
    integer, intent(in) :: itb
    integer             :: inp

    call system_clock(start)

    inp = size(Aexc)

    select case(itb)
    case (11)
      allocate(bX1exc(inp))
      bX1exc(:) = Aexc(1:inp)
    case (12)
      allocate(bX2exc(inp))
      bX2exc(:) = Aexc(1:inp)
    case (21)
      allocate(bY1exc(inp))
      bY1exc(:) = Aexc(1:inp)
    case (22)
      allocate(bY2exc(inp))
      bY2exc(:) = Aexc(1:inp)
    case (31)
      allocate(bZ1exc(inp))
      bZ1exc(:) = Aexc(1:inp)
    case (32)
      allocate(bZ2exc(inp))
      bZ2exc(:) = Aexc(1:inp)
    end select
    call system_clock(finish)
    call addTime(' bc', finish - start)
  end subroutine

  subroutine setBorderInside(itb, Aexc)
    integer, allocatable, intent(in) :: Aexc(:)
    integer, intent(in) :: itb
    integer             :: inp

    call system_clock(start)

    inp = size(Aexc)

    select case(itb)
    case (11)
      allocate(bX1ins(inp))
      bX1ins(:) = Aexc(1:inp)
    case (12)
      allocate(bX2ins(inp))
      bX2ins(:) = Aexc(1:inp)
    case (21)
      allocate(bY1ins(inp))
      bY1ins(:) = Aexc(1:inp)
    case (22)
      allocate(bY2ins(inp))
      bY2ins(:) = Aexc(1:inp)
    case (31)
      allocate(bZ1ins(inp))
      bZ1ins(:) = Aexc(1:inp)
    case (32)
      allocate(bZ2ins(inp))
      bZ2ins(:) = Aexc(1:inp)
    end select
    call system_clock(finish)
    call addTime(' bc', finish - start)
  end subroutine setBorderInside
end module BC
