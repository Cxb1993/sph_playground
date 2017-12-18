module timing

  use arrayresize, only: resize

  implicit none

  public :: printTimes, addTime, init, destroy

  private
  save
    integer, parameter :: namelen = 15

    integer(8), allocatable  :: timings(:)
    character(len=namelen), allocatable :: names(:)
    real :: clockrate

  contains

  subroutine init()
    integer(8) :: cr, cm
    external rtc

    call system_clock(count_rate=cr, count_max=cm)
    clockrate = real(cr)
  end subroutine

  subroutine destroy()
    deallocate(names)
    deallocate(timings)
  end subroutine

  subroutine addTime(inkey, inval)
    character(len=*), intent(in) :: inkey
    integer(8), intent(in)      :: inval

    integer(8) :: nonzeroval
    integer :: n, i, f

    if (inval == 0) then
      ! huge error could be here, as it is hardware dependent
      nonzeroval = 1000
    else
      nonzeroval = inval
    end if

    n = size(names)
    f = 0

    if ( allocated(names) ) then
      do i = 1,n
        if ( names(i) == inkey ) then
          timings(i) = timings(i) + nonzeroval
          f = 1
        end if
      end do
      if ( f == 0 ) then
        call resize(names, namelen, n, n+1)
        call resize(timings, n, n+1)
        names(n+1)   = inkey
        timings(n+1) = nonzeroval
      end if
    else
      allocate(names(1))
      allocate(timings(1))
      names(1)   = inkey
      timings(1) = nonzeroval
    end if
  end subroutine

  subroutine printTimes()
    integer :: i
    real    :: totaltime, thistime

    if (allocated(names)) then
      totaltime = 0.
      print *, '##############################################'
      print *, '#####    Times:'
      do i = 1,size(timings)
        thistime = timings(i)/clockrate
        write(*, "(A, A, A, F15.6)") " # # ", names(i), ": ", thistime
        totaltime = totaltime + thistime
      end do
      write(*, "(A, F15.6)") " # #  total time    : ", totaltime
      print *, '##############################################'
    end if
  end subroutine
end module
