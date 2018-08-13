module rhofuncs
 implicit none

contains
  function MTIHopkins2017(y) result(f)
    real, intent(in)  :: y
    real              :: f, g, h
    h = 3.
    g = -1.
    f = 2.*(h - y)**(1. - g*h)
  end function MTIHopkins2017

  function MTILowresHopkins2017(y) result(f)
    real, intent(in)  :: y
    real              :: f, g, h
    g = -0.75
    h = 1. - g*0.5
    ! f = g*y + h
    if (y < 0) f = -1./10.**6/1.5/5.*(-5. - y)**9
    if (y > 0) f = -1./10.**6/1.5/5.*(-5. + y)**9
  end function MTILowresHopkins2017

  function uniRho(y) result(f)
    real, intent(in)  :: y
    real :: f
    f = y
  end function uniRho
end module
