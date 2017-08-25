module kernel
 use physcon, only:pi
 implicit none
 character(len=11), public :: kernelname = 'zipM6'
 real, parameter, public  :: radkern  = 2
 real, parameter, public  :: radkern2 = 4
 real, parameter, public  :: cnormk = 0.0281249999999998/pi
 real, parameter, public  :: wab0 = 66., gradh0 = -3.*wab0
 real, parameter, public  :: dphidh0 = 1.44040178571444
 real, parameter, public  :: cnormk_drag = 0.0452008928571463/pi
 real, parameter, public  :: hfact_default = 1.3

contains

pure subroutine get_kernel(q2,q,wkern,grkern)
 real, intent(in)  :: q2,q
 real, intent(out) :: wkern,grkern
 real :: q4

 if (1.5*q < 1) then
    q4 = q2*q2
    wkern  = -75.9375*q4*q + 151.875*q4 - 135.*q2 + 66.
    grkern = q*(-379.6875*q2*q + 607.5*q2 - 270.)
 elseif (1.5*q < 2) then
    wkern  = -(1.5*q - 3)**5 + 6*(1.5*q - 2)**5
    grkern = -7.5*(1.5*q - 3)**4 + 45.*(1.5*q - 2)**4
 elseif (1.5*q < 3) then
    wkern  = -(1.5*q - 3)**5
    grkern = -7.5*(1.5*q - 3)**4
 else
    wkern  = 0
    grkern = 0
 endif

end subroutine

pure elemental real function wkern(q2,q)
 real, intent(in) :: q2,q
 real :: q4

 if (1.5*q < 1) then
    q4 = q2*q2
    wkern = -75.9375*q4*q + 151.875*q4 - 135.*q2 + 66.
 elseif (1.5*q < 2) then
    wkern = -(1.5*q - 3)**5 + 6*(1.5*q - 2)**5
 elseif (1.5*q < 3) then
    wkern = -(1.5*q - 3)**5
 else
    wkern = 0
 endif

end function

pure elemental real function grkern(q2,q)
 real, intent(in) :: q2,q

 if (1.5*q < 1) then
    grkern = q*(-379.6875*q2*q + 607.5*q2 - 270.)
 elseif (1.5*q < 2) then
    grkern = -7.5*(1.5*q - 3)**4 + 45.*(1.5*q - 2)**4
 elseif (1.5*q < 3) then
    grkern = -7.5*(1.5*q - 3)**4
 else
    grkern = 0
 endif

end function

pure subroutine get_kernel_grav1(q2,q,wkern,grkern,dphidh)
 real, intent(in)  :: q2,q
 real, intent(out) :: wkern,grkern,dphidh
 real :: q4, q6

 if (1.5*q < 1) then
    q4 = q2*q2
    q6 = q4*q2
    wkern  = -75.9375*q4*q + 151.875*q4 - 135.*q2 + 66.
    grkern = q*(-379.6875*q2*q + 607.5*q2 - 270.)
    dphidh = 1.22042410714285*q6*q - 2.84765624999998*q6 + 3.79687499999997*q4 - &
                 3.71249999999997*q2 + 1.44040178571444
 elseif (1.5*q < 2) then
    q4 = q2*q2
    q6 = q4*q2
    wkern  = -(1.5*q - 3)**5 + 6*(1.5*q - 2)**5
    grkern = -7.5*(1.5*q - 3)**4 + 45.*(1.5*q - 2)**4
    dphidh = -0.610212053571424*q6*q + 4.27148437499997*q6 - 11.3906249999999*q4*q + &
                 13.2890624999999*q4 - 4.21874999999997*q2*q - 2.86874999999998*q2 &
                 + 1.42533482142873
 elseif (1.5*q < 3) then
    q4 = q2*q2
    q6 = q4*q2
    wkern  = -(1.5*q - 3)**5
    grkern = -7.5*(1.5*q - 3)**4
    dphidh = 0.122042410714285*q6*q - 1.42382812499999*q6 + 6.83437499999995*q4*q - &
                 17.0859374999999*q4 + 22.7812499999998*q2*q - 13.6687499999999*q2 &
                 + 2.19676339285722
 else
    wkern  = 0
    grkern = 0
    dphidh = 0
 endif

end subroutine

pure subroutine kernel_softening(q2,q,potensoft,fsoft)
 real, intent(in)  :: q2,q
 real, intent(out) :: potensoft,fsoft
 real :: q4, q6

 if (1.5*q < 1) then
    q4 = q2*q2
    q6 = q4*q2
    potensoft = -0.152553013392856*q6*q + 0.406808035714282*q6 - 0.759374999999994*q4 + &
                 1.23749999999999*q2 - 1.44040178571444 - 0.894531249999644/q
    fsoft     = -1.06787109374999*q6 + 2.44084821428569*q4*q - 3.03749999999998*q2*q + &
                 2.47499999999998*q + 0.894531249999644/q2
 elseif (1.5*q < 2) then
    q4 = q2*q2
    q6 = q4*q2
    potensoft = 0.076276506696428*q6*q - 0.610212053571424*q6 + 1.89843749999998*q4*q - &
                 2.65781249999998*q4 + 1.05468749999999*q2*q + &
                 0.956249999999992*q2 - 1.42533482142873 - 0.894845145088926/q
    fsoft     = 0.533935546874996*q6 - 3.66127232142854*q4*q + 9.49218749999992*q4 - &
                 10.6312499999999*q2*q + 3.16406249999997*q2 + 1.91249999999998*q &
                 + 0.894845145088926/q2
 elseif (1.5*q < 3) then
    q4 = q2*q2
    q6 = q4*q2
    potensoft = -0.0152553013392856*q6*q + 0.203404017857141*q6 - 1.13906249999999*q4*q &
                 + 3.41718749999997*q4 - 5.69531249999995*q2*q + &
                 4.55624999999996*q2 - 2.19676339285722 - 0.862702287946196/q
    fsoft     = -0.106787109374999*q6 + 1.22042410714285*q4*q - 5.69531249999995*q4 + &
                 13.6687499999999*q2*q - 17.0859374999999*q2 + 9.11249999999993*q &
                 + 0.862702287946196/q2
 else
    potensoft = -1/q
    fsoft     = q**(-2)
 endif

end subroutine

!------------------------------------------
! double-humped version of the kernel for
! use in drag force calculations
!------------------------------------------
pure elemental real function wkern_drag(q2,q)
 real, intent(in) :: q2,q
 real :: q4

 if (1.5*q < 1) then
    q4 = q2*q2
    wkern_drag = q2*(-75.9375*q4*q + 151.875*q4 - 135.*q2 + 66.)
 elseif (1.5*q < 2) then
    wkern_drag = q2*(-(1.5*q - 3)**5 + 6*(1.5*q - 2)**5)
 elseif (1.5*q < 3) then
    wkern_drag = -q2*(1.5*q - 3)**5
 else
    wkern_drag = 0
 endif

end function

end module
