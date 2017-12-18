module BC
  use const
  use timing,           only: addTime

  implicit none

  public :: fixed1, fixed3, destroy, &
            setBorder, set_sqare_box_sides,&
            getSqaureBoxSides, setBorderInside, periodic3v2, periodic1v2,&
            periodicV3

  public :: particlesnumber, bordersize, &
            xmin, xmax, ymin, ymax, zmin, zmax

  private
  save
    integer              :: nx, ny, nz
    integer, allocatable :: bX1ins(:), bY1ins(:), bZ1ins(:), bX2ins(:), bY2ins(:), bZ2ins(:)
    integer, allocatable :: bX1exc(:), bY1exc(:), bZ1exc(:), bX2exc(:), bY2exc(:), bZ2exc(:)
    integer(8)           :: start=0, finish=0
    ! GLOABL =(
    integer :: particlesnumber, bordersize
    real :: xmin, xmax, ymin, ymax, zmin, zmax

contains
  subroutine destroy
    if (allocated(bX1exc)) then
      deallocate(bX1exc)
    end if
    if (allocated(bX2exc)) then
      deallocate(bX2exc)
    end if
    if (allocated(bY1exc)) then
      deallocate(bY1exc)
    end if
    if (allocated(bY2exc)) then
      deallocate(bY2exc)
    end if
    if (allocated(bZ1exc)) then
      deallocate(bZ1exc)
    end if
    if (allocated(bZ2exc)) then
      deallocate(bZ2exc)
    end if

    if (allocated(bX1ins)) then
      deallocate(bX1ins)
    end if
    if (allocated(bX2ins)) then
      deallocate(bX2ins)
    end if
    if (allocated(bY1ins)) then
      deallocate(bY1ins)
    end if
    if (allocated(bY2ins)) then
      deallocate(bY2ins)
    end if
    if (allocated(bZ1ins)) then
      deallocate(bZ1ins)
    end if
    if (allocated(bZ2ins)) then
      deallocate(bZ2ins)
    end if
  end subroutine

  subroutine set_sqare_box_sides(inx, iny, inz)
    integer, intent(in) :: inx, iny, inz
    nx = inx
    ny = iny
    nz = inz
  end subroutine set_sqare_box_sides

  subroutine getSqaureBoxSides(onx, ony, onz)
    integer, intent(out) :: onx, ony, onz
    onx = nx
    ony = ny
    onz = nz
  end subroutine getSqaureBoxSides

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

  subroutine periodicV3(axis, A1, A2, A3, AA1, AA2, AA3)
    real, allocatable, optional, intent(inout) :: A1(:), A2(:), A3(:), AA1(:,:), AA2(:,:), AA3(:,:)
    integer, intent(in) :: axis
    integer             :: i, sbx, sby, sbz, smax

    call system_clock(start)

    select case(axis)
    case (ebc_all)
      sbx = size(bX1ins)
      sby = size(bY1ins)
      sbz = size(bZ1ins)
      smax = max(sbx, sby, sbz)
      do i =1,smax
        if (i <= sbx) then
          A1(bX1exc(i)) = A1(bX1ins(i))
          A1(bX2exc(i)) = A1(bX2ins(i))
        end if
        if (i <= sby) then
          A1(bY1exc(i)) = A1(bY1ins(i))
          A1(bY2exc(i)) = A1(bY2ins(i))
        end if
        if (i <= sbz) then
          A1(bZ1exc(i)) = A1(bZ1ins(i))
          A1(bZ2exc(i)) = A1(bZ2ins(i))
        end if
      end do
    end select
    call system_clock(finish)
    call addTime(' bc', finish - start)
  end subroutine periodicV3

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
    real, intent(out)   :: A(particlesnumber)

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
    real, intent(out)   :: A(3,particlesnumber)

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
