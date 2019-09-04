module sts_integrator
implicit none

public sts_integrate
real, save, allocatable, dimension(:) :: &
  sts_Y0,sts_Yjm1,sts_Yjm2,sts_Yj,sts_Ys,sts_DY0,sts_DYjm1

contains

real function b(j)
  integer, intent(in) :: j
  if (j <= 2) then
    b = 1./3.
  else
    b = (j*j + j - 2.)/(2.*j*(j + 1.))
  endif
end function b

real function a(j)
  integer, intent(in) :: j
  a = 1 - b(j)
end function a

real function mu(j)
  integer, intent(in) :: j
  mu = (2.*j - 1.)/j*b(j)/b(j-1)
end function mu

real function nu(j)
  integer, intent(in) :: j
  nu = -(j - 1.)/j*b(j)/b(j-2)
end function nu

real function mu_hat(j,s)
  integer, intent(in) :: j,s
  if (j == 1) then
    mu_hat = 4./(3.*(s*s + s - 2.))
  else
    mu_hat = 4.*(2.*j - 1.)/(j*(s*s + s - 2.))*b(j)/b(j-1)
  endif
end function mu_hat

real function gm_hat(j,s)
  integer, intent(in) :: j,s
  gm_hat = - a(j-1)*mu_hat(j,s)
end function gm_hat

subroutine test()
  print*, '        mu(2) = ', mu(2), ' err = ', mu(2) - 3./2.
  print*, '        mu(3) = ', mu(3), ' err = ', mu(3) - 25./12.

  print*, '        nu(2) = ', nu(2), ' err = ', nu(2) - (-1./2.)
  print*, '        nu(3) = ', nu(3), ' err = ', nu(3) - (-5./6.)

  print*, '  mu_hat(1,3) = ', mu_hat(1,3), ' err = ',mu_hat(1,3) - 2./15.
  print*, '  mu_hat(2,3) = ', mu_hat(2,3), ' err = ',mu_hat(2,3) - 3./5.
  print*, '  mu_hat(3,3) = ', mu_hat(3,3), ' err = ',mu_hat(3,3) - 5./6.

  print*, '  gm_hat(2,3) = ', gm_hat(2,3), ' err = ',gm_hat(2,3) - (-2./5.)
  print*, '  gm_hat(3,3) = ', gm_hat(3,3), ' err = ',gm_hat(3,3) - (-5./9.)

  print*, '1-mu(2)-nu(2) = ', 1.-mu(2)-nu(2), ' err = ',1-mu(2)-nu(2) - (0.)
  print*, '1-mu(3)-nu(3) = ', 1.-mu(3)-nu(3), ' err = ',1-mu(3)-nu(3) - (-1./4.)
end subroutine test

! subroutine get_sts_status()
!
! end subroutine get_sts_status

subroutine init(n)
  integer, intent(in) :: n

  allocate(sts_Y0(n))
  allocate(sts_Yjm1(n))
  allocate(sts_Yjm2(n))
  allocate(sts_Yj(n))
  allocate(sts_Ys(n))
  allocate(sts_DY0(n))
  allocate(sts_DYjm1(n))
end subroutine init

subroutine sts_integrate(n,sts_s,store,dt)
  use const
  use iterator, only: iterate

  real, allocatable, intent(inout) :: store(:,:)
  integer, intent(in) :: n,sts_s
  real, intent(in)    :: dt

  integer :: sts_j
  real :: dedt

  sts_Y0(:) = store(es_t,:)
  call iterate(n, store, dedt)
  sts_DY0(:)  = store(es_ddt,:)
  sts_Yjm2(:) = sts_Y0(:)
  sts_Yjm1(:) = sts_Y0(:) + mu_hat(1,sts_s)*dt*sts_DY0(:)
  ! letzgered
  do sts_j = 2,sts_s
    store(es_t,:) = sts_Yjm1(:)
    call iterate(n, store, dedt)
    sts_DYjm1(:)  = store(es_ddt,:)

    sts_Yj(:)   = mu(sts_j)*sts_Yjm1(:) &
      + nu(sts_j)*sts_Yjm2(:) &
      + ( 1. - mu(sts_j) - nu(sts_j))*sts_Y0(:) &
      + mu_hat(sts_j,sts_s)*dt*sts_DYjm1(:) &
      + gm_hat(sts_j,sts_s)*dt*sts_DY0(:)
    ! next cicle preparation
    sts_Yjm2(:) = sts_Yjm1(:)
    sts_Yjm1(:) = sts_Yj(:)
    ! if (eqSet(eqs_radexch)==1) then
    !   call update_radenergy(n,store,dt)
    ! endif
  enddo
  store(es_t,:) = sts_Yj(:)
end subroutine sts_integrate

end module sts_integrator
