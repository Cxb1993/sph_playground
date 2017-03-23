module printer
  implicit none

  public :: Output, AppendLine, getTime

  private

  integer, save :: ifile = 0
  real :: start=0., finish=0., elapsed=0.


contains
  subroutine getTime(ot)
    real, intent(out) :: ot
    ot = elapsed
  end subroutine getTime

  subroutine Output(time, ptype, x, v, dv, m, den, slen, pres, ien, cf, err)
    real, allocatable, intent(in)    :: x(:,:), v(:,:), dv(:,:), m(:), err(:), &
                                        den(:), slen(:), pres(:), ien(:), cf(:)
    integer, allocatable, intent(in) :: ptype(:)
    real, intent(in)    :: time
    character (len=40)  :: fname
    integer :: iu, j, n

    call cpu_time(start)

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
    call cpu_time(finish)
    elapsed = elapsed + (finish - start)
  end subroutine Output

  subroutine AppendLine(n, A, fname)
    integer           :: n
    real, intent(in)  :: A(n)
    character (len=*), intent(in) :: fname
    integer :: iu

    call cpu_time(start)

    open(newunit=iu, file=fname, status='old', form='formatted', access='append')
    write(iu, *) A
    close(iu)
    call cpu_time(finish)
    elapsed = elapsed + (finish - start)
  end subroutine AppendLine
end module printer
