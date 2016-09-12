module printer
  implicit none
  integer, save :: ifile = 0

contains

  subroutine output(pos, vel, acc, mas, den, slen, time, p_num)
    integer, intent(in) :: p_num
    real, intent(in)  :: pos(0:p_num), mas(0:p_num), vel(0:p_num), den(0:p_num), slen(0:p_num), acc(0:p_num)
    real, intent(in) :: time
    character (len=40) :: fname
    integer :: iu, j

    write(fname, "(a,i2.2)") 'step_', ifile
    open(newunit=iu, file=fname, status='replace', form='formatted')
    write(iu,*) time
    do j = 0, p_num
      write(iu, *) pos(j), vel(j), acc(j), mas(j), den(j), slen(j)
    end do
    close(iu)
    ifile = ifile + 1
  end subroutine output

  subroutine plot_simple(x, y, fname)
    real, intent(in)  :: x, y
    character (len=*), intent(in) :: fname
    open(17, file=fname, status='replace', form='formatted')
    write(17, *) x, y
    close(17,status='keep')
  end subroutine plot_simple
end module printer
