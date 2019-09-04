module iterator
  use errprinter,       only: error
  use eos
  use const
  use circuit1
  use circuit2,         only: c2
  use BC,               only: reflectRealParticlesPeriodically, &
                              createPeriodicBorder, &
                              clearPeriodicParticles, &
                              createFixedBorders, &
                              findInsideBorderParticles, &
                              updateFixedToSymmetric
  use state,            only: getdim,&
                              get_equations, &
                              ginitvar,&
                              getPartNumber,&
                              getProcess,&
                              getGamma
  use neighboursearch,  only: findneighboursKDT

 implicit none

 public :: iterate

 private
 save
 real, allocatable  :: dbtmp(:,:)
 integer            :: initialised = 0, ivt, ipc

contains
  subroutine initIterate(n, store)
    real, allocatable, intent(inout) :: &
      store(:,:)
    integer, intent(in) :: &
      n

    allocate(dbtmp(3,n))
    call findneighboursKDT(store)
    call ginitvar(ivt)
    call getProcess(ipc)
    initialised = 1
  end subroutine

  subroutine iterate(n, store, maxconsenrg)
    real, allocatable, intent(inout) :: &
      store(:,:)
    integer, intent(in) :: n
    real, intent(out) :: maxconsenrg

    integer :: dim, ttp, rpn, fpn, i
    real    :: gamma

    call getdim(dim)
    call get_equations(ttp)
    call getPartNumber(r=rpn, f=fpn)
    call getGamma(gamma)

    if (initialised == 0) call initIterate(n, store)

    select case(ipc)
    case(epc_fullyperiodic)
      call clearPeriodicParticles(store)
      call reflectRealParticlesPeriodically(store, ebc_all)
      call findInsideBorderParticles(store)
      call createPeriodicBorder(store, ebc_all)

      call findneighboursKDT(store)
      call c1(store)
      call eos_adiabatic(store, gamma)
      call c2(store, maxconsenrg)
    case(epc_borderless)
      call findneighboursKDT(store)
      call c1(store)
      call eos_adiabatic(store, gamma)
      call c2(store, maxconsenrg)

    case(epc_backcompatibility)
      select case(ivt)
      case (ett_shock12)
        call clearPeriodicParticles(store)
        call findInsideBorderParticles(store)
        call createPeriodicBorder(store, ebc_y)

        call findneighboursKDT(store)
        call c1(store)
        call eos_adiabatic(store, gamma)
        call c2(store, maxconsenrg)
      case (ett_pulse, ett_ring, ett_fld_gauss)
        call findneighboursKDT(store)
        call c1(store)
        call c2(store, maxconsenrg)
      case (ett_soundwave)
        call clearPeriodicParticles(store)
        call reflectRealParticlesPeriodically(store, ebc_all)
        call findInsideBorderParticles(store)
        call createPeriodicBorder(store, ebc_all)
        call findneighboursKDT(store)
        call c1(store)
        ! call eos_adiabatic(store, gamma)
        call eos_isothermal(store)
        call c2(store, maxconsenrg)
      case (ett_hydroshock)
        call clearPeriodicParticles(store)
        call reflectRealParticlesPeriodically(store, ebc_y)
        call reflectRealParticlesPeriodically(store, ebc_z)
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
        call reflectRealParticlesPeriodically(store, ebc_x)
        call findInsideBorderParticles(store)
        call createPeriodicBorder(store, ebc_x)
        call findneighboursKDT(store)
        call c1(store)
        call eos_adiabatic(store, gamma)
        call c2(store, maxconsenrg)
        store(es_ay,1:rpn) = store(es_ay,1:rpn) - 1.
      case(ett_mtilowres)
        call clearPeriodicParticles(store)
        call reflectRealParticlesPeriodically(store, ebc_x)
        call findInsideBorderParticles(store)
        call createPeriodicBorder(store, ebc_x)
        call findneighboursKDT(store)
        call c1(store)
        call eos_adiabatic(store, gamma)
        call c2(store, maxconsenrg)
        do i = 1,rpn
          store(es_ay,i) = store(es_ay,i) - erf(store(es_ry,i) * 8.)
          ! print*, "if (store(es_ry,i) > 0.00001) store(es_ay,i) = store(es_ay,i) - 1."
          ! print*, "if (store(es_ry,i) < -0.00001) store(es_ay,i) = store(es_ay,i) + 1."
          ! print*, -erf(store(es_ry,i) * 8.)
          ! print*, store(es_ry,i)*8.
          ! print*, '----'
          ! read*
          ! if (store(es_ry,i) > 0.00001) store(es_ay,i) = store(es_ay,i) - 1.
          ! if (store(es_ry,i) < -0.00001) store(es_ay,i) = store(es_ay,i) + 1.
        end do
      case (ett_OTvortex)
        call clearPeriodicParticles(store)
        call reflectRealParticlesPeriodically(store, ebc_all)
        call findInsideBorderParticles(store)
        call createPeriodicBorder(store, ebc_all)
        call findneighboursKDT(store)
        call c1(store)
        call eos_adiabatic(store, gamma)
        call c2(store, maxconsenrg)
      case default
        call error('Task type was not set', '', __FILE__, __LINE__)
      end select
    case default
      call error('Process type was not set', '', __FILE__, __LINE__)
    end select
  end subroutine iterate
end module
