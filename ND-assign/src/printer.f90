module printer
  implicit none

  public :: print_output, print_appendline

  private

  integer, save :: ifile = 0

contains

  subroutine print_output(time, ptype, x, v, dv, m, den, slen, pres, ien, cf, err)
    real, allocatable, intent(in)    :: x(:,:), v(:,:), dv(:,:), m(:), err(:), &
                                        den(:), slen(:), pres(:), ien(:), cf(:)
    integer, allocatable, intent(in) :: ptype(:)
    real, intent(in)    :: time
    character (len=40)  :: fname
    integer :: iu, j, n
    n = size(ptype)
    write(fname, "(a,i5.5)") 'output/step_', ifile
    open(newunit=iu, file=fname, status='replace', form='formatted')
    write(iu,*) time
    do j = 1, n
      if (ptype(j) == 1) then
        write(iu, *) x(:,j), v(:,j), dv(:,j), m(j), den(j), slen(j), pres(j), ien(j), cf(j), err(j)
      end if
    end do
    close(iu)
    ifile = ifile + 1
  end subroutine print_output

  subroutine print_appendline(n, A, fname)
    integer           :: n
    real, intent(in)  :: A(n)
    character (len=*), intent(in) :: fname
    integer :: iu

    open(newunit=iu, file=fname, status='old', form='formatted', access='append')
    write(iu, *) A
    close(iu)
  end subroutine print_appendline
end module printer
