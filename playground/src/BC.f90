module BC
  implicit none

  public :: periodic1, periodic3, fixed1, fixed3
  public :: set_particles_numbers, set_border, set_sqare_box_sides, getSqaureBoxSides

  private
    integer, save              :: ns, nbrd, nx, ny, nz
    integer, allocatable, save :: borderX1(:), borderY1(:), borderZ1(:), borderX2(:), borderY2(:), borderZ2(:)

contains
  subroutine set_particles_numbers(ins, inbrd)
    integer, intent(in) :: ins, inbrd
    ns = ins
    nbrd = inbrd
  end subroutine set_particles_numbers

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

  subroutine set_border(itb, inp, A)
    integer, intent(in) :: itb, inp, A(inp)

    if(inp.ne.0) then
      select case(itb)
      case (11)
        allocate(borderX1(inp))
        borderX1 = A
      case (12)
        allocate(borderX2(inp))
        borderX2 = A
      case (21)
        allocate(borderY1(inp))
        borderY1 = A
      case (22)
        allocate(borderY2(inp))
        borderY2 = A
      case (31)
        allocate(borderZ1(inp))
        borderZ1 = A
      case (32)
        allocate(borderZ2(inp))
        borderZ2 = A
      end select
    end if
  end subroutine set_border
  !
  !-- See page 19
  !
  subroutine periodic1(A, axis)
    integer, intent(in) :: axis
    real, intent(out)   :: A(ns)
    integer             :: i

    select case(axis)
    case (1)
      ! A(borderX1) = A(borderX2 - (nbrd + 1) * ny * nz)
      ! A(borderX2) = A(borderX1 + (nbrd + 1) * ny * nz)
      do i = 1, size(borderX1)
        A(borderX1(i)) = A(borderX2(i) - (nbrd + 1) * ny * nz)
        A(borderX2(i)) = A(borderX1(i) + (nbrd + 1) * ny * nz)
      end do
    case (2)
      ! A(borderY1) = A(borderY2 - (nbrd + 1) * nz)
      ! A(borderY2) = A(borderY1 + (nbrd + 1) * nz)
      do i = 1, size(borderY1)
        A(borderY1(i)) = A(borderY2(i) - (nbrd + 1) * nz)
        A(borderY2(i)) = A(borderY1(i) + (nbrd + 1) * nz)
      end do
    case (3)
      ! A(borderZ1) = A(borderZ2 - (nbrd + 1))
      ! A(borderZ2) = A(borderZ1 + (nbrd + 1))
      do i = 1, size(borderZ1)
        A(borderZ1(i)) = A(borderZ2(i) - (nbrd + 1))
        A(borderZ2(i)) = A(borderZ1(i) + (nbrd + 1))
      end do
    end select
  end subroutine periodic1

  subroutine periodic3(A, axis, dim)
    integer, intent(in) :: axis, dim
    real, intent(out)   :: A(3,ns)
    integer             :: i

    select case(axis)
    case (00)
      if (size(borderX1) /= 1) then
        do i = 1, size(borderX1)
          A(:,borderX1(i)) = A(:,borderX2(i) - (nbrd + 1) * ny * nz)
          A(:,borderX2(i)) = A(:,borderX1(i) + (nbrd + 1) * ny * nz)
        end do
      end if
      if (dim > 1) then
        if (size(borderY1) /= 1) then
          do i = 1, size(borderY1)
            A(:,borderY1(i)) = A(:,borderY2(i) - (nbrd + 1) * nz)
            A(:,borderY2(i)) = A(:,borderY1(i) + (nbrd + 1) * nz)
          end do
        end if
      end if
      if (dim == 3) then
        if (size(borderZ1) /= 1) then
          do i = 1, size(borderZ1)
            A(:,borderZ1(i)) = A(:,borderZ2(i) - (nbrd + 1))
            A(:,borderZ2(i)) = A(:,borderZ1(i) + (nbrd + 1))
          end do
        end if
      end if
    case (10)
      do i = 1, size(borderX1)
        A(:,borderX1(i)) = A(:,borderX2(i) - (nbrd + 1) * ny * nz)
        A(:,borderX2(i)) = A(:,borderX1(i) + (nbrd + 1) * ny * nz)
      end do
    case (20)
      do i = 1, size(borderY1)
        A(:,borderY1(i)) = A(:,borderY2(i) - (nbrd + 1) * nz)
        A(:,borderY2(i)) = A(:,borderY1(i) + (nbrd + 1) * nz)
      end do
    case (30)
      do i = 1, size(borderZ1)
        A(:,borderZ1(i)) = A(:,borderZ2(i) - (nbrd + 1))
        A(:,borderZ2(i)) = A(:,borderZ1(i) + (nbrd + 1))
      end do
    end select
  end subroutine periodic3

  subroutine fixed1(A, axeside, k)
    integer, intent(in) :: axeside
    real, intent(in)    :: k
    real, intent(out)   :: A(ns)

    select case(axeside)
    case (00)
      A(borderX1) = k
      A(borderX2) = k
      A(borderY1) = k
      A(borderY2) = k
      A(borderZ1) = k
      A(borderZ2) = k
    case (11)
      A(borderX1) = k
    case (12)
      A(borderX2) = k
    case (21)
      A(borderY1) = k
    case (22)
      A(borderY2) = k
    case (31)
      A(borderZ1) = k
    case (32)
      A(borderZ2) = k
    end select
  end subroutine fixed1

  subroutine fixed3(A, axeside, dim, k)
    integer, intent(in) :: axeside, dim
    real, intent(in)    :: k
    real, intent(out)   :: A(3,ns)

    select case(axeside)
    case (00)
      A(dim,borderX1) = k
      A(dim,borderX2) = k
      A(dim,borderY1) = k
      A(dim,borderY2) = k
      A(dim,borderZ1) = k
      A(dim,borderZ2) = k
    case (11)
      A(dim,borderX1) = k
    case (12)
      A(dim,borderX2) = k
    case (21)
      A(dim,borderY1) = k
    case (22)
      A(dim,borderY2) = k
    case (31)
      A(dim,borderZ1) = k
    case (32)
      A(dim,borderZ2) = k
    end select
  end subroutine fixed3
end module BC
