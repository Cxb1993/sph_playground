module printer
  implicit none

  public :: print_output, plot_simple

  private

  integer, save :: ifile = 0

contains

  subroutine print_output(n, time, pos, vel, acc, mas, den, slen, pres, ien, cf, err)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n,3), mas(n), vel(n,3), acc(n,3),&
                           den(n), slen(n), pres(n), ien(n), cf(n), err(n)
    real, intent(in)    :: time
    character (len=40)  :: fname
    integer :: iu, j

    write(fname, "(a,i5.5)") 'output/step_', ifile
    open(newunit=iu, file=fname, status='replace', form='formatted')
    write(iu,*) time
    do j = 1, n
      write(iu, *) pos(j,:), vel(j,:), acc(j,:), mas(j), den(j), slen(j), pres(j), ien(j), cf(j), err(j)
    end do
    close(iu)
    ifile = ifile + 1
  end subroutine print_output

  subroutine plot_simple(x, n, A, fname)
    integer           :: n
    real, intent(in)  :: x, A(n)
    character (len=*), intent(in) :: fname
    integer :: iu

    open(newunit=iu, file=fname, status='old', form='formatted', access='append')
    write(iu, *) x, A
    close(iu)
  end subroutine plot_simple
end module printer
