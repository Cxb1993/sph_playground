module rad_exchange
 implicit none
 public :: update_radenergy!,set_radfluxesandregions

 private

contains

subroutine update_radenergy(n,store,dt)
  use const
  use errprinter,       only: warning, error
  use state,            only: getStateVal

  integer, intent(in) ::n
  real, intent(in)    ::dt
  real, intent(inout) ::store(:,:)

  real :: ui,pmassi,rhoi,xii
  real :: ack,a,cv1,kappa,dudt,etot,unew,gamma,gmw
  integer :: i

  call getStateVal(ec_gamma,gamma)
  call getStateVal(ec_molecularmass,gmw)

  a   = 4.*sigmaB/lightspeed
  cv1 = (gamma-1.)*gmw/Rg

  !$omp parallel do default(none)&
  !$omp private(kappa,ack,rhoi,ui)&
  !$omp private(dudt,xii,etot,unew)&
  !$omp shared(store)&
  !$omp shared(pmassi)&
  !$omp shared(n,dt,cv1,a)
  do i = 1,n
    pmassi = store(es_m,i)
    kappa  = store(es_kappa,i)
    ack    = 4.*sigmaB*kappa

    rhoi = store(es_den,i)
    ui   = store(es_u,i)
    dudt = store(es_du,i)
    xii  = store(es_t,i)
    etot = ui + xii
    unew = ui
    if (xii <= 0) then
       call warning('radiation energy is negative before exchange',i,__FILE__,__LINE__)
    endif
    ! if (i==3273) then
    !    print*, 'Before:  ', 'T_gas=',unew*cv1,'T_rad=',(rhoi*(etot-unew)/a)**(1./4.)
    ! endif
    call solve_internal_energy_implicit(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    ! call solve_internal_energy_implicit_substeps(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt)
    ! call solve_internal_energy_explicit_substeps(unew,ui,rhoi,etot,dudt,ack,a,cv1,dt,di)
    store(es_u,i) = unew
    store(es_t,i) = etot - unew
    ! if (i==3273) then
    !     print*, 'After:   ', 'T_gas=',unew*cv1,'T_rad=',(rhoi*(etot-unew)/a)**(1./4.)
    !     read*
    ! endif
    if (store(es_t,i) <= 0) &
       call error('radiation energy is negative after exchange',i,__FILE__,__LINE__)
    if (store(es_u,i) <= 0) &
       call error('thermal energy is negative after exchange',i,__FILE__,__LINE__)
  end do
  !$omp end parallel do
end subroutine update_radenergy

subroutine solve_internal_energy_implicit_substeps(unew,ui,rho,etot,dudt,ack,a,cv1,dt)
 real, intent(out) :: unew
 real, intent(in)  :: ui, etot, dudt, dt, rho, ack, a, cv1

 real     :: fu,dfu,eps,dts,dunew,unewp,uip
 integer  :: iter,i,level

 unew = ui
 unewp = 2.*ui
 eps = 1e-8
 dunew = (unewp-unew)/unew
 level = 1
 do while((abs(dunew) > eps).and.(level <= 2**10))
   unewp = unew
   dts   = dt/level
   uip   = ui
   do i=1,level
     iter  = 0
     fu    = huge(1.)
     do while((abs(fu) > eps).and.(iter < 10))
       iter = iter + 1
       fu   = unew/dts - uip/dts - dudt - ack*(rho*(etot-unew)/a - (unew*cv1)**4)
       dfu  = 1./dts + ack*(rho/a + 4.*(unew**3*cv1**4))
       unew = unew - fu/dfu
     end do
     uip = unew
   end do
   dunew = (unewp-unew)/unew
   level = level*2
 end do
end subroutine solve_internal_energy_implicit_substeps

subroutine solve_internal_energy_implicit(unew,u0,rho,etot,dudt,ack,a,cv1,dt)
 real, intent(out) :: unew
 real, intent(in)  :: u0, etot, dudt, dt, rho, ack, a, cv1
 real     :: fu,dfu,eps,uold
 integer  :: iter

 unew = u0
 uold = 2*u0

 iter = 0
 eps = 1e-16
 do while ((abs(unew-uold) > eps).and.(iter < 10))
   uold = unew
   iter = iter + 1
   fu   = unew/dt - u0/dt - dudt - ack*(rho*(etot-unew)/a - (unew*cv1)**4)
   dfu  = 1./dt + ack*(rho/a + 4.*(unew**3*cv1**4))
   unew = unew - fu/dfu
 end do
end subroutine solve_internal_energy_implicit

subroutine solve_internal_energy_explicit(unew,ui,rho,etot,dudt,ack,a,cv1,dt)
 real, intent(out) :: unew
 real, intent(in)  :: ui, etot, dudt, dt, rho, ack, a, cv1

 unew = ui + dt*(dudt + ack*(rho*(etot-ui)/a - (ui*cv1)**4))
end subroutine solve_internal_energy_explicit

subroutine solve_internal_energy_explicit_substeps(unew,ui,rho,etot,dudt,ack,a,cv1,dt)
 real, intent(out) :: unew
 real, intent(in)  :: ui, etot, dudt, dt, rho, ack, a, cv1
 real     :: du,eps,dts,unews,uis
 integer  :: i,level

 level = 1
 eps   = 1e-8
 dts   = dt

 unew = ui + dts*(dudt + ack*(rho*(etot-ui)/a - (ui*cv1)**4))
 unews = 2*unew
 du = (unew-unews)/unew
 do while((abs(du) > eps).and.(level <= 2**10))
   unews = unew
   level = level*2
   dts  = dt/level
   unew = 0.
   uis  = ui
   do i=1,level
     unew = uis + dts*(dudt + ack*(rho*(etot-uis)/a - (uis*cv1)**4))
     uis  = unew
   end do
   du = (unew-unews)/unew
 end do
end subroutine solve_internal_energy_explicit_substeps

! subroutine set_radfluxesandregions(npart,radiation,xyzh,vxyzu)
!   use part,    only: igas,massoftype,rhoh,ifluxx,ifluxy,ifluxz,ithick,iradxi,ikappa
!   use eos,     only: get_spsound
!   use options, only: ieos
!   use physcon, only:c
!   use units,   only:unit_velocity
!
!   real, intent(inout)    :: radiation(:,:),vxyzu(:,:)
!   real, intent(in)       :: xyzh(:,:)
!   integer, intent(in)    :: npart
!
!   integer :: i
!   real :: pmassi,H,cs,rhoi,r,c_code,lambdai,prevfrac
!
!   pmassi = massoftype(igas)
!
!   prevfrac = 100.*count(radiation(ithick,:)==1)/real(size(radiation(ithick,:)))
!   sch = sch * (1. + 0.005 * (80. - prevfrac))
!   print*, "-}+{- RADIATION Scale Height Multiplier: ", sch
!
!   radiation(ithick,:) = 1
!
!   c_code = c/unit_velocity
!
!   do i = 1,npart
!     rhoi = rhoh(xyzh(4,i),pmassi)
!     ! if (rhoi < 2e-4) then
!     !   if (xyzh(1,i) < 0.) then
!     !     radiation(ifluxy:ifluxz,i) = 0.
!     !     radiation(ifluxx,i)        = -rhoi*abs(radiation(iradxi,i))*0.5
!     !     radiation(ithick,i) = 0
!     !   else if (xyzh(1,i) > 0.) then
!     !     radiation(ifluxy:ifluxz,i) = 0.
!     !     radiation(ifluxx,i)        =  rhoi*abs(radiation(iradxi,i))*0.5
!     !     radiation(ithick,i) = 0
!     !   end if
!     ! end if
!     cs = get_spsound(ieos,xyzh(:,i),rhoi,vxyzu(:,i))
!     r  = sqrt(dot_product(xyzh(1:3,i),xyzh(1:3,i)))
!     H  = cs*sqrt(r**3)
!     if (abs(xyzh(3,i)) > sch*H) then
!        radiation(ithick,i) = 0
!        ! if (xyzh(3,i) <= 0.) then
!        !   radiation(ifluxx:ifluxy,i) = 0.
!        !   radiation(ifluxz,i)        =  rhoi*abs(radiation(iradxi,i))
!        !   radiation(ithick,i) = 0
!        ! else if (xyzh(3,i) > 0.) then
!        !   radiation(ifluxx:ifluxy,i) = 0.
!        !   radiation(ifluxz,i)        =  -rhoi*abs(radiation(iradxi,i))
!        !   radiation(ithick,i) = 0
!        ! end if
!     end if
!     ! Ri = sqrt(&
!     !   dot_product(radiation(ifluxx:ifluxz,i),radiation(ifluxx:ifluxz,i)))&
!     !   /(radiation(ikappa,i)*rhoi*rhoi*radiation(iradxi,i))
!     ! lambda = (2. + Ri)/(6. + 3*Ri + Ri*Ri)
!     ! lambdai = 1./3
!     ! ! print*, xyzh(4,i), lambdai/radiation(ikappa,i)/rhoi*c_code/cs
!     ! ! read*
!     ! if (xyzh(4,i) > lambdai/radiation(ikappa,i)/rhoi*c_code/cs) then
!     !    radiation(ithick,i) = 0
!     ! endif
!   end do
! end subroutine set_radfluxesandregions
! subroutine mcfost_do_analysis()
!
! end subroutine mcfost_do_analysis
!
end module rad_exchange
