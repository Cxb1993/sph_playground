module iterator
  use eos
  use const
  use circuit1
  use circuit2,         only: c2
  use BC,               only: reflecParticlesPeriodic, &
                              createPhantomPeriodic, &
                              clearPhantomParticles, &
                              createPhantomFixed, &
                              realpartnumb, &
                              artpartnumb, &
                              needcrosref
  use state,            only: get_difftype,&
                              getdim,&
                              get_tasktype, &
                              ginitvar
  use neighboursearch,  only: findneighboursKDT

 implicit none

 public :: iterate

 private
 save
 real, allocatable  :: dbtmp(:,:)
 integer            :: initialised = 0

contains
  subroutine initIterate(n)
    integer, intent(in) :: n

    allocate(dbtmp(3,n))
    initialised = 1
  end subroutine

  subroutine iterate(n, hfac, gamma, store)
    real, allocatable, intent(inout) :: &
      store(:,:)
    integer, intent(in) :: n
    real, intent(in)    :: hfac, gamma

    integer             :: dim, ttp, dtp, ivt

    call getdim(dim)
    call get_tasktype(ttp)
    call get_difftype(dtp)
    call ginitvar(ivt)

    if (initialised == 0) then
      call initIterate(n)
    end if

    select case(ivt)
    case (ett_pulse, ett_ring)
      ! call clearPhantomParticles(store)
      call createPhantomPeriodic(store, ebc_all)
      call findneighboursKDT(store)
      call c1(store, hfac)
      call c2(store)
    ! case (ett_soundwave)
    !   call clearPhantomParticles(pos)
    !   call reflecParticlesPeriodic(pos, ebc_all)
    !   call createPhantomPeriodic(pos, ebc_all)
    !   call findneighboursKDT_V2(pos, h)
    !   call c1(pos, mas, sk, h, den, om, cf, dcf)
    !   ! call eos_adiabatic(den, uei, prs, c, gamma)
    !   call eos_isothermal(den, c(1), prs)
    !   call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
    ! case (ett_hydroshock)
    !   call clearPhantomParticles(pos)
    !   call createPhantomFixed(pos, ebc_x)
    !   call createPhantomPeriodic(pos, ebc_y)
    !   call createPhantomPeriodic(pos, ebc_z)
    !   call findneighboursKDT_V2(pos, h)
    !   call c1(pos, mas, sk, h, den, om, cf, dcf)
    !   call eos_adiabatic(den, uei, prs, c, gamma)
    !   call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
    ! case (ett_alfvenwave)
    !   call findneighboursKDT(ptype, pos, h)
    !
    !   call c1(pos, mas, sk, h, den, om, cf, dcf)
    !   call periodic1v2(den, ebc_all)
    !   call periodic1v2(h,   ebc_all)
    !   call periodic1v2(om,  ebc_all)
    !
    !   call eos_adiabatic(den, uei, prs, c, gamma)
    !   ! call eos_isothermal(den, c(1), prs)
    !   call periodic1v2(prs, ebc_all)
    !   call periodic1v2(c,   ebc_all)
    !
    !   call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
    !   call periodic3v2(acc, ebc_all)
    !   call periodic3v2(dcf, ebc_all)
    !   call periodic1v2(due, ebc_all)
    !   call periodic1v2(dh,  ebc_all)
    !   dbtmp(:,:) = kcf(:,2,:)
    !   call periodic3v2(dbtmp, ebc_all)
    !   kcf(:,2,:) = dbtmp(:,:)
    ! case(ett_mti)
    !   call clearPhantomParticles(pos)
    !   call reflecParticlesPeriodic(pos, ebc_x)
    !   call createPhantomPeriodic(pos, ebc_x)
    !   artpn1 = artpartnumb
    !   call createPhantomFixed(pos, ebc_y)
    !   artpn2 = artpartnumb
    !   ! expand all storage and fill [artpn1:artpn2] with what ever I want
    !   print*, size(pos,dim=2), realpartnumb, artpn1, artpn2
    !   read*
    !   call findneighboursKDT_V2(pos, h)
    !   call c1(pos, mas, sk, h, den, om, cf, dcf)
    !   call eos_adiabatic(den, uei, prs, c, gamma)
    !   call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
    ! case (ett_OTvortex)
    !   call clearPhantomParticles(pos)
    !   call reflecParticlesPeriodic(pos, ebc_all)
    !   call createPhantomPeriodic(pos, ebc_all)
    !   call findneighboursKDT_V2(pos, h)
    !   call c1(pos, mas, sk, h, den, om, cf, dcf)
    !   call eos_adiabatic(den, uei, prs, c, gamma)
    !   call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
    case default
      print *, 'Task type was not defined in iterator:145.'
      stop
    end select
  end subroutine iterate
end module

! select case (ttp)
! case (1, 2, 3, 4, 9)
!   ! mooved to ivt check
! case(5)
!   ! 'diff-laplace'
!   print*, "FIX ME. I should depend on IVT not EQS"
!   call findneighboursN2plusStatic(ptype, pos, h)
!   select case(dtp)
!   case(1)
!     call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
!   case(2)
!     call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
!     call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
!   case default
!     print *, 'Diff type is not set in iterator'
!     stop
!   end select
! case(6)
!   ! 'diff-graddiv'
!   print*, "FIX ME. I should depend on IVT not EQS"
!   call findneighboursN2plusStatic(ptype, pos, h)
!   call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
!   call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
! case(7,8)
! case(10)
!   ! diff-artvisc
!   print*, "FIX ME. I should depend on IVT not EQS"
!   ! call findneighboursN2plus(ptype, pos, h)
!   call findneighboursKDT(ptype, pos, h)
!   call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
! case default
!   print *, 'Task type was not defined in iterator.f90: line 140.'
!   stop
! end select
