module printer
  implicit none
  integer, save :: ifile = 0

contains

  subroutine output(n, time, pos, vel, acc, mas, den, slen, pres, ien)
    integer, intent(in) :: n
    real, intent(in)    :: pos(3,n), mas(n), vel(n), den(n), slen(n), acc(n), pres(n), ien(n)
    real, intent(in)    :: time
    character (len=40)  :: fname
    integer :: iu, j

    write(fname, "(a,i5.5)") 'steps/output_', ifile
    open(newunit=iu, file=fname, status='replace', form='formatted')
    write(iu,*) time
    do j = 1, n
      write(iu, *) pos(1,j), pos(2,j), pos(3,j), vel(j), acc(j), mas(j), den(j), slen(j), pres(j), ien(j)
    end do
    close(iu)
    ifile = ifile + 1
  end subroutine output

  subroutine plot_simple(x, y, fname)
    real, intent(in)  :: x, y
    character (len=*), intent(in) :: fname
    integer :: iu

    open(newunit=iu, file=fname, status='old', form='formatted', access='append')
    write(iu, *) x, y
    close(iu)
  end subroutine plot_simple
end module printer
