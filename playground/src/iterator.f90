module iterator
  use eos
  use const
  use circuit1
  use circuit2,         only: c2
  use BC,               only: reflecPeriodicParticles, &
                              createPeriodicBorder, &
                              clearPeriodicParticles, &
                              createFixedBorders, &
                              findInsideBorderParticles, &
                              updateFixedToSymmetric
  use state,            only: getdim,&
                              get_equations, &
                              ginitvar,&
                              getPartNumber
  use neighboursearch,  only: findneighboursKDT, &
                              findneighboursN2

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

  subroutine iterate(n, gamma, store, maxconsenrg)
    real, allocatable, intent(inout) :: &
      store(:,:)
    integer, intent(in) :: n
    real, intent(in)    :: gamma
    real, intent(out) :: maxconsenrg

    integer             :: dim, ttp, ivt, rpn, fpn, i

    call getdim(dim)
    call get_equations(ttp)
    call ginitvar(ivt)
    call getPartNumber(r=rpn, f=fpn)

    if (initialised == 0) then
      call initIterate(n)
    end if

    select case(ivt)
    case (ett_shock12)
      call clearPeriodicParticles(store)
      call findInsideBorderParticles(store)
      call createPeriodicBorder(store, ebc_y)

      call findneighboursKDT(store)
      call c1(store)
      call c2(store, maxconsenrg)
    case (ett_pulse, ett_ring)
      call findneighboursKDT(store)
      call c1(store)
      call c2(store, maxconsenrg)
    case (ett_soundwave)
      call clearPeriodicParticles(store)
      call reflecPeriodicParticles(store, ebc_all)
      call findInsideBorderParticles(store)
      call createPeriodicBorder(store, ebc_all)
      call findneighboursKDT(store)
      call c1(store)
      ! call eos_adiabatic(store, gamma)
      call eos_isothermal(store)
      call c2(store, maxconsenrg)
    case (ett_hydroshock)
      call clearPeriodicParticles(store)
      call reflecPeriodicParticles(store, ebc_y)
      call reflecPeriodicParticles(store, ebc_z)
      call findInsideBorderParticles(store)
      call createPeriodicBorder(store, ebc_y)
      call createPeriodicBorder(store, ebc_z)
      call findneighboursKDT(store)
      call c1(store)
      call eos_adiabatic(store, gamma)
      call c2(store, maxconsenrg)
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
    case(ett_mti)
      call clearPeriodicParticles(store)
      call reflecPeriodicParticles(store, ebc_x)
      call findInsideBorderParticles(store)
      call createPeriodicBorder(store, ebc_x)
      call findneighboursKDT(store)
      call c1(store)
      call eos_adiabatic(store, gamma)
      call c2(store, maxconsenrg)
      store(es_ay,1:rpn) = store(es_ay,1:rpn) - 1.
    case(ett_mtilowres)
      call clearPeriodicParticles(store)
      call reflecPeriodicParticles(store, ebc_x)
      call findInsideBorderParticles(store)
      call createPeriodicBorder(store, ebc_x)
      call findneighboursKDT(store)
      call c1(store)
      call eos_adiabatic(store, gamma)
      call c2(store, maxconsenrg)
      do i = 1,rpn
        if (store(es_ry,i) > 0) store(es_ay,i) = store(es_ay,i) - 1.
        if (store(es_ry,i) < 0) store(es_ay,i) = store(es_ay,i) + 1.
      end do
    case (ett_OTvortex)
      call clearPeriodicParticles(store)
      call reflecPeriodicParticles(store, ebc_all)
      call findInsideBorderParticles(store)
      call createPeriodicBorder(store, ebc_all)
      call findneighboursKDT(store)
      call c1(store)
      call eos_adiabatic(store, gamma)
      call c2(store, maxconsenrg)
    case default
      print *, 'Task type was not defined in ./src/iterator.f90:140'
      stop
    end select
  end subroutine iterate
end module
