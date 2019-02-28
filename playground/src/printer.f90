module printer
  use const

  use kernel, only: get_krad
  use timing, only: addTime
  use state,  only: getPartNumber,&
                    getStateVal,&
                    setStateVal

  implicit none

  public :: outputAscii, AppendLine, outputScreen, handleOutput

  private
  integer(8), save :: start=0, finish=0

contains

  subroutine handleOutput(iter, n, t, sumdt, dedt, dt, sumdedt, store, err)
    use dumper, only: dump

    integer, intent(in)::iter, n
    real, intent(in) :: t, sumdt, dedt, dt, sumdedt
    real, allocatable, intent(in) :: store(:,:), err(:)

    integer :: lastnpic, silent, usedumps

    call getStateVal(ec_silent, silent)
    call getStateVal(ec_usedumps, usedumps)

    if (silent == 0) call outputAscii(t, store, err)
    if (usedumps == 1) call dump(store, t)
    call outputScreen(iter, n, t, sumdt, dedt, dt, sumdedt, store)
    if ((silent==0).or.(usedumps==1)) then
      call getStateVal(ec_lastprint, lastnpic)
      call setStateVal(ec_lastprint, real(lastnpic+1))
    end if
  end subroutine handleOutput

  subroutine outputScreen(iter, n, t, sumdt, dedt, dt, sumdedt, store)
    integer, intent(in)::iter, n
    real, intent(in) :: t, sumdt, dedt, dt, sumdedt
    real, allocatable, intent(in) :: store(:,:)

    integer :: lastprintnumber

    call getStateVal(ec_lastprint, lastprintnumber)

    write(*, fmt="(A, I7, A, ES7.1, A, ES10.4, A, ES10.4, A, ES10.4, A, ES10.4, A, A, ES10.4, A, ES10.4)") &
      " #", lastprintnumber, &
      " | i=", real(iter), &
      " | t=", t, &
      " | dt=", sumdt/(iter+1), &
      " | h=[", minval(store(es_h,1:n)), ":", maxval(store(es_h,1:n)), "]",&
      " | dedt=", dedt*dt,&
      " | S(dedt)=", sumdedt
  end subroutine outputScreen

  subroutine outputAscii(time, store, err)
    real, allocatable, intent(in) :: &
      store(:,:), err(:)
    real, intent(in)    :: time
    real                :: kr, e
    character (len=40)  :: fname
    integer :: iu = 0, j, n, rn, ifile

    call system_clock(start)
    call get_krad(kr)
    call getPartNumber(r=rn)
    call getStateVal(ec_lastprint, ifile)

    n = size(store,2)
    write(fname, "(a,i5.5)") 'output/step_', ifile
    open(newunit=iu, file=fname, status='replace', form='formatted')
    write(iu,*) time
    do j = 1, n
      e = 0.
      if (j <= rn) e = err(j)
      ! if (int(store(es_type,j)) /= ept_empty) then
      if (int(store(es_type,j)) /= ept_empty) then
      ! if (int(store(es_type,j)) == ept_real) then
        write(iu, *) store(es_rx:es_rz,j), store(es_vx:es_vz,j),&
          store(es_ax:es_az,j),store(es_m,j),store(es_den,j),&
          store(es_h,j),store(es_p,j),store(es_u,j),store(es_t,j),&
          store(es_dtdx:es_dtdz,j), e
      end if
    end do
    close(iu)
    call system_clock(finish)
    call addTime(' printer', finish - start)
    ! print*, 1
  end subroutine outputAscii

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
