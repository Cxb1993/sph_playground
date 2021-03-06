!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2017 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: table_utils
!
!  DESCRIPTION:
!   This module contains low-level utilities related to tabulated data
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id: cb93c13f923f89e3d60b216c87a09a1a26be8f22 $
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module tables
 implicit none

 public :: yinterp, linspace

 private

contains

!----------------------------------------------------
!+
!  Function to linearly interpolate values from
!  from an array
!+
!----------------------------------------------------
pure real function yinterp(ytab,xtab,x)
 real, intent(in) :: ytab(:),xtab(:)
 real, intent(in) :: x
 integer :: i,n
 real    :: xprev,xi,dx,dydx

 n = size(xtab)

 yinterp = 0.
 if (n <= 0) return

 if (x <= xtab(1)) then
    yinterp = ytab(1)
    return
 elseif (x >= xtab(n)) then
    yinterp = ytab(n)
    return
 endif

 xprev = xtab(1)
 over_tab: do i=2,n
    xi = xtab(i)
    if (x >= xprev .and. x < xi) then
       dx      = x - xprev
       dydx    = (ytab(i) - ytab(i-1))/(xi - xprev)
       yinterp = ytab(i-1) + dx*dydx
       exit over_tab
    endif
    xprev = xi
 enddo over_tab

end function yinterp

!--------------------------------------------------------
!+
!  Function to fill an array with equally spaced points
!+
!--------------------------------------------------------
pure subroutine linspace(x,xmin,xmax)
 real, intent(out) :: x(:)
 real, intent(in)  :: xmin,xmax
 integer :: i, n
 real    :: dx

 n  = size(x)
 dx = (xmax - xmin)/real(n-1)
 do i=1,n
    x(i) = xmin + (i-1)*dx
 enddo

end subroutine linspace

end module tables
