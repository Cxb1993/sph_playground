module timing

  use utils, only: resize

  implicit none

  public :: printTimes, addTime, init

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

  call system_clock(count_rate=cr)
  call system_clock(count_max=cm)
  clockrate = real(cr)

  end subroutine init

  subroutine addTime(inkey, inval)
    character(len=*), intent(in) :: inkey
    integer(8), intent(in)      :: inval

    integer :: n, i, f

    n = size(names)
    f = 0

    if ( allocated(names) ) then
      do i = 1,n
        if ( names(i) == inkey ) then
          timings(i) = timings(i) + inval
          f = 1
        end if
      end do
      if ( f == 0 ) then
        call resize(names, namelen, n, n+1)
        call resize(timings, n, n+1)
        names(n+1)   = inkey
        timings(n+1) = inval
      end if
    else
      allocate(names(1))
      allocate(timings(1))
      names(1)   = inkey
      timings(1) = inval
    end if
  end subroutine addTime

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
  end subroutine printTimes
end module timing
