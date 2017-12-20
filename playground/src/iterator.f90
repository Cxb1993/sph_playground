module iterator
  use eos
  use const
  use circuit1
  use circuit2,         only: c2
  use BC,               only: periodic1v2,&
                              periodic3v2,&
                              reflecParticlesPeriodic, &
                              createPhantomPeriodic
  use state,            only: get_difftype,&
                              getdim,&
                              get_tasktype, &
                              ginitvar
  use neighboursearch,  only: findneighboursN2plusStatic, &
                              findneighboursKDT,&
                              findneighboursKDT_V2

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

  subroutine iterate(n, sk, gamma, ptype, pos, vel, acc, &
                    mas, den, h, dh, om, prs, c, uei, due, cf, dcf, kcf)
    real, allocatable, intent(inout), dimension(:,:)  :: pos, vel, acc, cf, dcf
    real, allocatable, intent(inout), dimension(:)    :: mas, den, dh, prs, c, uei, due, om, h
    real, allocatable, intent(inout), dimension(:,:,:):: kcf
    integer, allocatable, intent(in)   :: ptype(:)

    integer, intent(in) :: n
    real, intent(in)    :: sk, gamma
    integer             :: dim, ttp, dtp, ivt, i

    call getdim(dim)
    call get_tasktype(ttp)
    call get_difftype(dtp)
    call ginitvar(ivt)

    dcf(:,:) = 0.
    if (initialised == 0) then
      call initIterate(n)
    end if

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
    select case(ivt)
    case (ett_ring)
      call findneighboursKDT(ptype, pos, h)

      call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
      call periodic1v2(den, ebc_all)
      call periodic1v2(h,   ebc_all)
      call periodic1v2(om,  ebc_all)
      ! call periodic3v2(dcf, ebc_all)

      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      ! call periodic3v2(dcf, ebc_all)
      call periodic1v2(dh,  ebc_all)
      ! do I need to do it
      kcf(:,2,:) = 0.
      acc(:,:) = 0.
    case (ett_soundwave)
      call findneighboursKDT(ptype, pos, h)

      call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
      call periodic1v2(den, ebc_all)
      call periodic1v2(h,   ebc_all)
      call periodic1v2(om,  ebc_all)

      ! call eos_adiabatic(den, uei, prs, c, gamma)
      call eos_isothermal(den, c(1), prs)
      call periodic1v2(prs, ebc_all)
      call periodic1v2(c,   ebc_all)

      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      call periodic3v2(acc, ebc_all)
      call periodic3v2(dcf, ebc_all)
      call periodic1v2(due, ebc_all)
      call periodic1v2(dh,  ebc_all)
    case (ett_hydroshock)
      call findneighboursKDT(ptype, pos, h)
      call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
      if ( dim > 1 ) then
        call periodic1v2(den, ebc_y)
        call periodic1v2(h,   ebc_y)
        call periodic1v2(om,  ebc_y)
      end if
      call eos_adiabatic(den, uei, prs, c, gamma)
      if ( dim > 1 ) then
        call periodic1v2(prs, ebc_y)
        call periodic1v2(c,   ebc_y)
      end if
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      if ( dim > 1 ) then
        call periodic3v2(acc, ebc_y)
        call periodic1v2(due, ebc_y)
        call periodic1v2(dh,  ebc_y)
      end if
    case (ett_alfvenwave)
      call findneighboursKDT(ptype, pos, h)

      call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
      call periodic1v2(den, ebc_all)
      call periodic1v2(h,   ebc_all)
      call periodic1v2(om,  ebc_all)

      call eos_adiabatic(den, uei, prs, c, gamma)
      ! call eos_isothermal(den, c(1), prs)
      call periodic1v2(prs, ebc_all)
      call periodic1v2(c,   ebc_all)

      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      call periodic3v2(acc, ebc_all)
      call periodic3v2(dcf, ebc_all)
      call periodic1v2(due, ebc_all)
      call periodic1v2(dh,  ebc_all)
      dbtmp(:,:) = kcf(:,2,:)
      call periodic3v2(dbtmp, ebc_all)
      kcf(:,2,:) = dbtmp(:,:)
    case(ett_mti)
      call findneighboursKDT(ptype, pos, h)

      call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)

      call periodic1v2(den, ebc_x)
      call periodic1v2(h,   ebc_x)
      call periodic1v2(om,  ebc_x)

      call eos_adiabatic(den, uei, prs, c, gamma)
      ! call eos_isothermal(den, c(1), prs)
      call periodic1v2(prs, ebc_x)
      call periodic1v2(c,   ebc_x)

      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      call periodic3v2(acc, ebc_x)
      call periodic3v2(dcf, ebc_x)
      call periodic1v2(due, ebc_x)
      call periodic1v2(dh,  ebc_x)
      ! call fixed3(acc, ebc_y, ebc_all, 0.)
      ! call fixed3(dcf, ebc_y, ebc_all, 0.)
      ! call fixed1(due, ebc_y, 0.)
      ! call fixed1(dh,  ebc_y, 0.)
      dbtmp(:,:) = kcf(:,2,:)
      call periodic3v2(dbtmp, ebc_x)
      ! call fixed3(dbtmp, ebc_y, ebc_all, 0.)
      kcf(:,2,:) = dbtmp(:,:)
    case (ett_OTvortex)
      call reflecParticlesPeriodic(pos)
      call createPhantomPeriodic(pos)
      call findneighboursKDT_V2(ptype, pos, h)
      call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
      call eos_adiabatic(den, uei, prs, c, gamma)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
    case default
      print *, 'Task type was not defined in iterator.f90: line 204.'
      stop
    end select
  end subroutine iterate
end module
