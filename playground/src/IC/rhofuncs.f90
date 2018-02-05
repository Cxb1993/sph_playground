module rhofuncs
 implicit none

 public :: MTIHopkins2017

contains
  function MTIHopkins2017(y) result(f)
    real, intent(in)  :: y
    real              :: f, g, h
    h = 3.
    g = -1.
    f = (h - y)**(1. - g*h)
    ! f = 1.
    ! f = (0.76 + (0.1-y)*1.28*10)
  end function MTIHopkins2017
end module
