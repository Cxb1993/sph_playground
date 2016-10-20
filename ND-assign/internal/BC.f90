module BC
  use kernel

  implicit none

  public :: set_periodic, set_fixed1, set_fixed3, set_periodic_dim
  public :: set_ns, set_borders, set_tasktype, get_tasktype

  private
    character (len=40)          :: tasktype
    integer, save              :: ns
    integer, allocatable, save :: borderX(:), borderY(:), borderZ(:)

contains
  subroutine set_ns(ins)
    integer, intent(in) :: ins
    ns = ins
  end subroutine set_ns

  subroutine set_borders(nx, ny, nz, ibx, iby, ibz)
    integer, intent(in) :: nx, ny, nz, ibx(nx), iby(ny), ibz(nz)

    allocate(borderX(nx))
    allocate(borderY(ny))
    allocate(borderZ(nz))

    borderX(:) = ibx(:)
    borderY(:) = iby(:)
    borderZ(:) = ibz(:)
  end subroutine set_borders

  subroutine set_tasktype(itt)
    character (len=*), intent(in) :: itt
    tasktype = itt
  end subroutine set_tasktype

  subroutine get_tasktype(ott)
    character (len=*), intent(out) :: ott
    ott = tasktype
  end subroutine get_tasktype

  subroutine set_periodic(n, A)
    integer, intent(in) :: n
    real, intent(out)   :: A(n)
    integer             :: i, nr, bn
    bn = 2

    nr = n - 2 * bn
    do i = 1, bn
      A(i) = A(nr + i)
      A(nr + bn + i) = A(bn + i)
    end do
  end subroutine set_periodic

  subroutine set_fixed1(A)
    real, intent(out) :: A(ns)
    integer           :: dim

    call get_dim(dim)

    A(borderX) = 0.
    if(dim.gt.1) then
      A(borderY) = 0.
      if(dim.eq.3) then
        A(borderZ) = 0.
      end if
    end if
  end subroutine set_fixed1

  subroutine set_fixed3(A)
    real, intent(out) :: A(ns,3)
    integer           :: dim

    call get_dim(dim)

    A(borderX,1) = 0.
    if(dim.gt.1) then
      A(borderY,2) = 0.
      if(dim.eq.3) then
        A(borderZ,3) = 0.
      end if
    end if
  end subroutine set_fixed3

  subroutine set_periodic_dim(A, k)
    integer, intent(in) :: k
    real, intent(out)   :: A(ns,3)
    integer             :: dim

    call get_dim(dim)

    select case(k)
    case (1)
      A(borderX,1) = 0.
    case (2)
      A(borderY,2) = 0.
    case (3)
      A(borderZ,3) = 0.
    end select
  end subroutine set_periodic_dim

end module BC
