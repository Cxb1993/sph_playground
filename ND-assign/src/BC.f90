module BC
  implicit none

  public :: set_periodic1, set_periodic3, set_fixed1, set_fixed3
  public :: set_particles_numbers, set_border

  private
    integer, save              :: ns, nbrd
    integer, allocatable, save :: borderX1(:), borderY1(:), borderZ1(:), borderX2(:), borderY2(:), borderZ2(:)

contains
  subroutine set_particles_numbers(ins, inbrd)
    integer, intent(in) :: ins, inbrd
    ns = ins
    nbrd = inbrd
  end subroutine set_particles_numbers

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

  subroutine set_periodic1(A, axe)
    integer, intent(in) :: axe
    real, intent(out)   :: A(ns)

    select case(axe)
    case (1)
      A(borderX1) = A(borderX2 - nbrd)
      A(nbrd + borderX1) = A(borderX2)
    case (2)
      A(borderY1) = A(borderY2 - nbrd)
      A(nbrd + borderY1) = A(borderY2)
    case (3)
      A(borderZ1) = A(borderZ2 - nbrd)
      A(nbrd + borderZ1) = A(borderZ2)
    end select
  end subroutine set_periodic1

  subroutine set_periodic3(A, axe, dim)
    integer, intent(in) :: axe, dim
    real, intent(out)   :: A(ns,3)

    select case(axe)
    case (10)
      A(borderX1,:) = A(borderX2 - nbrd,:)
      A(nbrd + borderX1,:) = A(borderX2,:)
    case (1)
      A(borderX1,dim) = A(borderX2 - nbrd,dim)
      A(nbrd + borderX1,dim) = A(borderX2,dim)
    case (20)
      A(borderY1,:) = A(borderY2 - nbrd,:)
      A(nbrd + borderY1,:) = A(borderY2,:)
    case (2)
      A(borderY1,dim) = A(borderY2 - nbrd,dim)
      A(nbrd + borderY1,dim) = A(borderY2,dim)
    case (3)
      A(borderZ1,dim) = A(borderZ2 - nbrd,dim)
      A(nbrd + borderZ1,dim) = A(borderZ2,dim)
    end select
  end subroutine set_periodic3

  subroutine set_fixed1(A, axeside, k)
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
  end subroutine set_fixed1

  subroutine set_fixed3(A, axeside, dim, k)
    integer, intent(in) :: axeside, dim
    real, intent(in)    :: k
    real, intent(out)   :: A(ns,3)

    select case(axeside)
    case (00)
      A(borderX1, dim) = k
      A(borderX2, dim) = k
      A(borderY1, dim) = k
      A(borderY2, dim) = k
      A(borderZ1, dim) = k
      A(borderZ2, dim) = k
    case (11)
      A(borderX1, dim) = k
    case (12)
      A(borderX2, dim) = k
    case (21)
      A(borderY1, dim) = k
    case (22)
      A(borderY2, dim) = k
    case (31)
      A(borderZ1, dim) = k
    case (32)
      A(borderZ2, dim) = k
    end select
  end subroutine set_fixed3

end module BC
