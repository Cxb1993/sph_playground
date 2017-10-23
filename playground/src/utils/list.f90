module list
  implicit none

  type intlist
    private
    integer, allocatable  :: elements(:)
    integer               :: length = 0
    integer               :: size = 0
  contains
    procedure :: append, toarr, print, e, xe, llen
  end type

contains

  subroutine append(list, element)
    use arrayresize, only: resize

    class(intlist)      :: list
    integer, intent(in) :: element

    if (list%size == 0) then
      allocate(list%elements(2))
      list%size = 2
      list%length = 1
      list%elements(1) = element
    else
      if (list%size == list%length) then
        call resize(list%elements, list%size, list%size*2)
        list%size = list%size*2
      end if
      list%length = list%length + 1
      list%elements(list%length) = element
    end if
  end subroutine

  pure function toarr(list) result(elements)
    class(intlist), intent(in)  :: list
    integer, allocatable        :: elements(:)

    allocate(elements(list%length))
    elements(:) = list%elements(1:list%length)
  end function

  subroutine print(list)
    class(intlist) :: list
    print*, ' Len: ', list%length, ' Size: ', list%size, ' List: ', list%elements(1:list%length)
  end subroutine

  pure subroutine e(list, idx, element)
    class(intlist), intent(in)  :: list
    integer, intent(in)         :: idx
    integer, intent(out)        :: element

    if ((idx > 0).and.(idx <= list%length)) then
      element = list%elements(idx)
    else
      error stop 'index is out of range in list'
    end if
  end subroutine

  pure subroutine xe(list, idx, element)
    class(intlist), intent(in)  :: list
    integer, intent(in)         :: idx
    integer, intent(out)        :: element

    element = list%elements(idx)
  end subroutine

  pure function llen(list) result(length)
    class(intlist), intent(in)  :: list
    integer                     :: length

    length = list%length
  end function
end module
