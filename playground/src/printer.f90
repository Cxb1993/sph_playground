module printer
  use const
  use kernel, only: get_krad
  use timing, only: addTime
  use state,  only: getPartNumber,&
                    getLastPrint,&
                    setLastPrint

  implicit none

  public :: Output, AppendLine

  private

  integer(8) :: start=0, finish=0


contains
  subroutine Output(time, store, err)
    real, allocatable, intent(in) :: &
      store(:,:), err(:)
    real, intent(in)    :: time
    real                :: kr, e
    character (len=40)  :: fname
    integer :: iu = 0, j, n, rn, ifile

    call system_clock(start)
    call get_krad(kr)
    call getPartNumber(r=rn)
    call getLastPrint(ifile)

    n = size(store,2)
    write(fname, "(a,i5.5)") 'output/step_', ifile
    open(newunit=iu, file=fname, status='replace', form='formatted')
    write(iu,*) time
    do j = 1, n
      e = 0.
      if (j <= rn) e = err(j)
      if (int(store(es_type,j)) /= ept_empty) then
      ! if (int(store(es_type,j)) /= ept_empty) then
      ! if (int(store(es_type,j)) == ept_real) then
        write(iu, *) store(es_rx:es_rz,j), store(es_vx:es_vz,j),&
          store(es_ax:es_az,j),store(es_m,j),store(es_den,j),&
          store(es_h,j),store(es_p,j),store(es_u,j),store(es_t,j),&
          store(es_bx:es_bz,j), e
      end if
    end do
    close(iu)
    call setLastPrint(ifile + 1)
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
