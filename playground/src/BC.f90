module BC
  use list
  use const
  use timing,           only: addTime

  implicit none

  public :: fixed1, fixed3, destroy, getCrossRef, &
            setBorder, setBorderInside, periodic3v2, periodic1v2, &
            reflecParticlesPeriodic, createPhantomPeriodic

  public :: realpartnumb, bordersize, artpartnumb, &
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
    integer :: realpartnumb=0, artpartnumb=0, initdone=0
    real :: xmin, xmax, ymin, ymax, zmin, zmax, bordersize

contains
  subroutine clearFastAllLists()
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
  end subroutine clearFastAllLists

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

  subroutine init()
    if (initdone == 0) then
      allocate(bX1ins(1))
      allocate(bY1ins(1))
      allocate(bZ1ins(1))
      allocate(bX2ins(1))
      allocate(bY2ins(1))
      allocate(bZ2ins(1))

      allocate(bX1exc(1))
      allocate(bY1exc(1))
      allocate(bZ1exc(1))
      allocate(bX2exc(1))
      allocate(bY2exc(1))
      allocate(bZ2exc(1))
      initdone = 1
    end if
  end subroutine init

  pure function getCrossRef(ini) result(oti)
    integer, intent(in) :: ini
    integer oti
    oti = posref%xe(ini)
  end function getCrossRef

  subroutine reflecParticlesPeriodic(pos)
    real, allocatable, intent(inout) :: pos(:,:)
    integer :: i
    real :: boxmax(3), boxmin(3), dxmin(3), dxmax(3)

    boxmax = [xmax, ymax, zmax]
    boxmin = [xmin, ymin, zmin]
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
  end subroutine reflecParticlesPeriodic

  subroutine createPhantomPeriodic(pos)
    use arrayresize,  only: resize
    use state,        only: getdim

    real, allocatable, intent(inout) :: pos(:,:)
    integer :: i, k, dim, spos, refidx
    real :: boxmax(3), boxmin(3)

    call clearFastAllLists()
    call getdim(dim)
    boxmax = [xmax, ymax, zmax] - bordersize
    boxmin = [xmin, ymin, zmin] + bordersize
    if (dim == 1) then
      do i=1,realpartnumb
        call posref%append(i)
        if (pos(1,i) < boxmin(1)) then
          call ibx1%append(i)
        else if (pos(1,i) > boxmax(1)) then
          call ibx2%append(i)
        end if
      end do
    else if (dim == 2) then
      do i=1,realpartnumb
        call posref%append(i)
        if (pos(1,i) < boxmin(1)) then
          call ibx1%append(i)
        else if (pos(1,i) > boxmax(1)) then
          call ibx2%append(i)
        end if
        if (pos(2,i) < boxmin(2)) then
          call iby1%append(i)
        else if (pos(2,i) > boxmax(2)) then
          call iby2%append(i)
        end if
      end do
    else if (dim == 3) then
      do i=1,realpartnumb
        call posref%append(i)
        if (pos(1,i) < boxmin(1)) then
          call ibx1%append(i)
        else if (pos(1,i) > boxmax(1)) then
          call ibx2%append(i)
        end if
        if (pos(2,i) < boxmin(2)) then
          call iby1%append(i)
        else if (pos(2,i) > boxmax(2)) then
          call iby2%append(i)
        end if
        if (pos(3,i) < boxmin(3)) then
          call ibz1%append(i)
        else if (pos(3,i) > boxmax(3)) then
          call ibz2%append(i)
        end if
      end do
    end if

    spos = size(pos,dim=2)
    artpartnumb = 1
    do i=1,ibx1%llen()
      k = realpartnumb+artpartnumb
      if (spos < k) then
        print*, "  <!> pos array was expanded due to phantom periodic particles. L:182"
        call resize(pos, k-1, 2*k)
        spos = 2*k
      end if
      call ebx1%append(k)
      call posref%append(ibx1%xe(i))
      pos(:,k) = pos(:,ibx1%xe(i))
      pos(1,k) = xmax + (pos(1,k) - xmin)
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
        print*, "  <!> pos array was expanded due to phantom periodic particles. L:208"
        call resize(pos, k-1, 2*k)
        spos = 2*k
      end if
      call ebx2%append(k)
      call posref%append(ibx2%xe(i))
      pos(:,k) = pos(:,ibx2%xe(i))
      pos(1,k) = xmin + (pos(1,k) - xmax)
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
    if (dim > 1) then
      do i=1,iby1%llen()
        k = realpartnumb+artpartnumb
        if (spos < k) then
          print*, "  <!> pos array was expanded due to phantom periodic particles. L:235"
          call resize(pos, k-1, 2*k)
          spos = 2*k
        end if
        call eby1%append(k)
        ! artificial particle can reference artificial again
        ! so I need to calculate dereference in while it is not a real particle
        refidx = iby1%xe(i)
        do while (refidx > realpartnumb)
          refidx = getCrossRef(refidx)
        end do
        call posref%append(refidx)
        pos(:,k) = pos(:,iby1%xe(i))
        pos(2,k) = ymax + (pos(2,k) - ymin)
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
          print*, "  <!> pos array was expanded due to phantom periodic particles. L:260"
          call resize(pos, k-1, 2*k)
          spos = 2*k
        end if
        call eby2%append(k)
        ! artificial particle can reference artificial again
        ! so I need to calculate dereference in while it is not a real particle
        refidx = iby2%xe(i)
        do while (refidx > realpartnumb)
          refidx = getCrossRef(refidx)
        end do
        call posref%append(refidx)
        pos(:,k) = pos(:,iby2%xe(i))
        pos(2,k) = ymin + (pos(2,k) - ymax)
        if (dim == 3) then
          if (pos(3,k) < boxmin(3)) then
            call ibz1%append(k)
          else if (pos(3,k) > boxmax(3)) then
            call ibz2%append(k)
          end if
        end if
        artpartnumb = artpartnumb + 1
      end do
      if (dim == 3) then
        do i=1,ibz1%llen()
          k = realpartnumb+artpartnumb
          if (spos < k) then
            print*, "  <!> pos array was expanded due to phantom periodic particles. L:286"
            call resize(pos, k-1, 2*k)
            spos = 2*k
          end if
          call ebz1%append(k)
          ! artificial particle can reference artificial again
          ! so I need to calculate dereference in while it is not a real particle
          refidx = ibz1%xe(i)
          do while (refidx > realpartnumb)
            refidx = getCrossRef(refidx)
          end do
          call posref%append(refidx)
          pos(:,k) = pos(:,ibz1%xe(i))
          pos(3,k) = zmax + (pos(3,k) - zmin)
          artpartnumb = artpartnumb + 1
        end do
        do i=1,ibz2%llen()
          k = realpartnumb+artpartnumb
          if (spos < k) then
            print*, "  <!> pos array was expanded due to phantom periodic particles. L:304"
            call resize(pos, k-1, 2*k)
            spos = 2*k
          end if
          call ebz2%append(k)
          ! artificial particle can reference artificial again
          ! so I need to calculate dereference in while it is not a real particle
          refidx = ibz2%xe(i)
          do while (refidx > realpartnumb)
            refidx = getCrossRef(refidx)
          end do
          call posref%append(refidx)
          pos(:,k) = pos(:,ibz2%xe(i))
          pos(3,k) = zmin + (pos(3,k) - zmax)
          artpartnumb = artpartnumb + 1
        end do
      end if
    end if
    artpartnumb = artpartnumb - 1
    ! call posref%print()
    ! print *, artpartnumb, realpartnumb, size(pos, dim=2)
    ! read*
  end subroutine createPhantomPeriodic

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
  end subroutine

  subroutine periodic1v2(A, axis)
    integer, intent(in)              :: axis
    real, allocatable, intent(inout) :: A(:)
    integer                          :: i

    call system_clock(start)

    select case(axis)
    case (ebc_all)
      do i = 1, size(bX1ins)
        A(bX1exc(i)) = A(bX1ins(i))
      end do
      do i = 1, size(bX2ins)
        A(bX2exc(i)) = A(bX2ins(i))
      end do
      do i = 1, size(bY1ins)
        A(bY1exc(i)) = A(bY1ins(i))
      end do
      do i = 1, size(bY2ins)
        A(bY2exc(i)) = A(bY2ins(i))
      end do
      do i = 1, size(bZ1ins)
        A(bZ1exc(i)) = A(bZ1ins(i))
      end do
      do i = 1, size(bZ2ins)
        A(bZ2exc(i)) = A(bZ2ins(i))
      end do
    case (ebc_x)
      do i = 1, size(bX1ins)
        A(bX1exc(i)) = A(bX1ins(i))
      end do
      do i = 1, size(bX2ins)
        A(bX2exc(i)) = A(bX2ins(i))
      end do
    case (ebc_y)
      do i = 1, size(bY1ins)
        A(bY1exc(i)) = A(bY1ins(i))
      end do
      do i = 1, size(bY2ins)
        A(bY2exc(i)) = A(bY2ins(i))
      end do
    case (ebc_z)
      do i = 1, size(bZ1ins)
        A(bZ1exc(i)) = A(bZ1ins(i))
      end do
      do i = 1, size(bZ2ins)
        A(bZ2exc(i)) = A(bZ2ins(i))
      end do
    end select
    call system_clock(finish)
    call addTime(' bc', finish - start)
  end subroutine periodic1v2

  ! subroutine periodicV3(axis, A1, A2, A3, AA1, AA2, AA3)
  !   real, allocatable, optional, intent(inout) :: A1(:), A2(:), A3(:), AA1(:,:), AA2(:,:), AA3(:,:)
  !   integer, intent(in) :: axis
  !   integer             :: i, sbx, sby, sbz, smax
  !
  !   call system_clock(start)
  !
  !   select case(axis)
  !   case (ebc_all)
  !     sbx = size(bX1ins)
  !     sby = size(bY1ins)
  !     sbz = size(bZ1ins)
  !     smax = max(sbx, sby, sbz)
  !     do i =1,smax
  !       if (i <= sbx) then
  !         A1(bX1exc(i)) = A1(bX1ins(i))
  !         A1(bX2exc(i)) = A1(bX2ins(i))
  !       end if
  !       if (i <= sby) then
  !         A1(bY1exc(i)) = A1(bY1ins(i))
  !         A1(bY2exc(i)) = A1(bY2ins(i))
  !       end if
  !       if (i <= sbz) then
  !         A1(bZ1exc(i)) = A1(bZ1ins(i))
  !         A1(bZ2exc(i)) = A1(bZ2ins(i))
  !       end if
  !     end do
  !   end select
  !   call system_clock(finish)
  !   call addTime(' bc', finish - start)
  ! end subroutine periodicV3

  subroutine periodic3v2(A, axis)
    integer, intent(in)              :: axis
    real, allocatable, intent(inout) :: A(:,:)
    integer                          :: i

    call system_clock(start)

    select case(axis)
    case (ebc_all)
      do i = 1, size(bX1ins)
        A(:,bX1exc(i)) = A(:,bX1ins(i))
      end do
      do i = 1, size(bX2ins)
        A(:,bX2exc(i)) = A(:,bX2ins(i))
      end do
      do i = 1, size(bY1ins)
        A(:,bY1exc(i)) = A(:,bY1ins(i))
      end do
      do i = 1, size(bY2ins)
        A(:,bY2exc(i)) = A(:,bY2ins(i))
      end do
      do i = 1, size(bZ1ins)
        A(:,bZ1exc(i)) = A(:,bZ1ins(i))
      end do
      do i = 1, size(bZ2ins)
        A(:,bZ2exc(i)) = A(:,bZ2ins(i))
      end do
    case (ebc_x)
      do i = 1, size(bX1ins)
        A(:,bX1exc(i)) = A(:,bX1ins(i))
      end do
      do i = 1, size(bX2ins)
        A(:,bX2exc(i)) = A(:,bX2ins(i))
      end do
    case (ebc_y)
      do i = 1, size(bY1ins)
        A(:,bY1exc(i)) = A(:,bY1ins(i))
      end do
      do i = 1, size(bY2ins)
        A(:,bY2exc(i)) = A(:,bY2ins(i))
      end do
    case (ebc_z)
      do i = 1, size(bZ1ins)
        A(:,bZ1exc(i)) = A(:,bZ1ins(i))
      end do
      do i = 1, size(bZ2ins)
        A(:,bZ2exc(i)) = A(:,bZ2ins(i))
      end do
    end select
    call system_clock(finish)
    call addTime(' bc', finish - start)
  end subroutine periodic3v2

  subroutine fixed1(A, axeside, k)
    integer, intent(in) :: axeside
    real, intent(in)    :: k
    real, intent(out)   :: A(realpartnumb)

    call system_clock(start)

    select case(axeside)
    case (ebc_all)
      A(bX1exc) = k
      A(bX2exc) = k
      A(bY1exc) = k
      A(bY2exc) = k
      A(bZ1exc) = k
      A(bZ2exc) = k
    case (ebc_x1)
      A(bX1exc) = k
    case (ebc_x2)
      A(bX2exc) = k
    case (ebc_x)
      A(bX1exc) = k
      A(bX2exc) = k
    case (ebc_y1)
      A(bY1exc) = k
    case (ebc_y2)
      A(bY2exc) = k
    case (ebc_y)
      A(bY1exc) = k
      A(bY2exc) = k
    case (ebc_z1)
      A(bZ1exc) = k
    case (ebc_z2)
      A(bZ2exc) = k
    case (ebc_z)
      A(bZ1exc) = k
      A(bZ2exc) = k
    end select
    call system_clock(finish)
    call addTime(' bc', finish - start)
  end subroutine fixed1

  subroutine fixed3(A, axeside, dim, k)
    integer, intent(in) :: axeside, dim
    real, intent(in)    :: k
    real, intent(out)   :: A(3,realpartnumb)

    call system_clock(start)

    select case(axeside)
    case (ebc_all)
      if (dim == ebc_all) then
        A(:,bX1exc) = k
        A(:,bX2exc) = k
        A(:,bY1exc) = k
        A(:,bY2exc) = k
        A(:,bZ1exc) = k
        A(:,bZ2exc) = k
      end if
      A(dim,bX1exc) = k
      A(dim,bX2exc) = k
      A(dim,bY1exc) = k
      A(dim,bY2exc) = k
      A(dim,bZ1exc) = k
      A(dim,bZ2exc) = k
    case (ebc_x1)
      A(dim,bX1exc) = k
      if (dim == ebc_all) then
        A(:,bX1exc) = k
      end if
    case (ebc_x2)
      A(dim,bX2exc) = k
      if (dim == ebc_all) then
        A(:,bX2exc) = k
      end if
    case (ebc_x)
      A(dim,bX1exc) = k
      A(dim,bX2exc) = k
      if (dim == ebc_all) then
        A(:,bX1exc) = k
        A(:,bX2exc) = k
      end if
    case (ebc_y1)
      A(dim,bY1exc) = k
      if (dim == ebc_all) then
        A(:,bY1exc) = k
      end if
    case (ebc_y2)
      A(dim,bY2exc) = k
      if (dim == ebc_all) then
        A(:,bY2exc) = k
      end if
    case (ebc_y)
      A(dim,bY1exc) = k
      A(dim,bY2exc) = k
      if (dim == ebc_all) then
        A(:,bY1exc) = k
        A(:,bY2exc) = k
      end if
    case (ebc_z1)
      A(dim,bZ1exc) = k
      if (dim == ebc_all) then
        A(:,bZ1exc) = k
      end if
    case (ebc_z2)
      if (dim == ebc_all) then
        A(:,bZ2exc) = k
      else
        A(dim,bZ2exc) = k
      end if
    case (ebc_z)
      if (dim == ebc_all) then
        A(:,bZ1exc) = k
        A(:,bZ2exc) = k
      else
        A(dim,bZ1exc) = k
        A(dim,bZ2exc) = k
      end if
    end select
    call system_clock(finish)
    call addTime(' bc', finish - start)
  end subroutine fixed3

  ! subroutine ReflectBorders(axis, pos)
  !   real, allocatable, optional, intent(inout) :: pos
  !   integer, intent(in) :: axis
  !   integer             :: i, sbx, sby, sbz, smax
  !
  !   call system_clock(start)
  !   select case(axis)
  !   case (ebc_all)
  !     sbx = size(bX1ins)
  !     sby = size(bY1ins)
  !     sbz = size(bZ1ins)
  !     smax = max(sbx, sby, sbz)
  !     do i =1,smax
  !       if (i <= sbx) then
  !         A1(bX1exc(i)) = A1(bX1ins(i))
  !         A1(bX2exc(i)) = A1(bX2ins(i))
  !       end if
  !       if (i <= sby) then
  !         A1(bY1exc(i)) = A1(bY1ins(i))
  !         A1(bY2exc(i)) = A1(bY2ins(i))
  !       end if
  !       if (i <= sbz) then
  !         A1(bZ1exc(i)) = A1(bZ1ins(i))
  !         A1(bZ2exc(i)) = A1(bZ2ins(i))
  !       end if
  !     end do
  !   end select
  !   call system_clock(finish)
  !   call addTime(' bc', finish - start)
  ! end subroutine ReflectBorders
end module BC
