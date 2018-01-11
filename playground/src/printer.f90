module printer
  use const
  use kernel, only: get_krad
  use BC,     only: realpartnumb,&
                    artpartnumb
  use timing, only: addTime

  implicit none

  public :: Output, AppendLine

  private

  integer, save :: ifile = 0
  integer(8) :: start=0, finish=0


contains
  subroutine Output(time, store, err)
    real, allocatable, intent(in) :: &
      store(:,:), err(:)
    real, intent(in)    :: time
    real                :: kr
    character (len=40)  :: fname
    integer :: iu = 0, j, n

    call system_clock(start)
    call get_krad(kr)


    n = realpartnumb + artpartnumb
    write(fname, "(a,i5.5)") 'output/step_', ifile
    open(newunit=iu, file=fname, status='replace', form='formatted')
    write(iu,*) time
    do j = 1, n
        write(iu, *) store(es_rx:es_rz,j), store(es_vx:es_vz,j),&
          store(es_ax:es_az,j),store(es_m,j),store(es_den,j),&
          store(es_h,j),store(es_p,j),store(es_u,j),store(es_t,j),&
          store(es_bx:es_bz,j), err(j)
    end do
    close(iu)
    ifile = ifile + 1
    call system_clock(finish)
    call addTime(' printer', finish - start)
  end subroutine Output

  subroutine AppendLine(A, fname, t)
    real, allocatable, intent(inout) :: A(:)
    character (len=*), intent(in) :: fname
    integer(8), intent(out) :: t
    integer :: iu = 0
    logical :: exist

    call system_clock(start)

    inquire(file=fname, exist=exist)
    if (exist) then
      open(newunit=iu, file=fname, status='old', form='formatted', access='append')
    else
      open(newunit=iu, file=fname, status='new', form='formatted')
    end if
    write(iu, *) A(:)
    close(iu)

    call system_clock(finish)
    t = finish - start
    call addTime(' printer', t)
  end subroutine AppendLine
end module printer
