module rhofuncs
 implicit none

 public :: MTIHopkins2017

contains
  pure function MTIHopkins2017(r) result(f)
    real, intent(in)  :: r
    real              :: f
    f = (1. - r/3.)**(-(1. + 1./3.)/(1./3.))
  end function MTIHopkins2017
end module
