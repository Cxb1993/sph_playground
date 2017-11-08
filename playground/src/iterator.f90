module iterator
  use eos
  use circuit1
  use timing,           only: addTime
  use circuit2,         only: c2, c15
  use BC
  use state,            only: get_difftype,&
                              getdim,&
                              get_tasktype, &
                              ginitvar
  use neighboursearch,  only: findneighboursN2plusStatic, &
                              findneighboursKDT

 implicit none

 public :: iterate

 private
 save
 integer(8)  :: start=0, finish=0

contains
  subroutine iterate(n, sk, gamma, ptype, pos, vel, acc, &
                    mas, den, h, dh, om, prs, c, uei, due, cf, dcf, kcf)
    real, allocatable, intent(inout), dimension(:,:)  :: pos, vel, acc, cf, dcf
    real, allocatable, intent(inout), dimension(:)    :: mas, den, dh, prs, c, uei, due, om, h
    real, allocatable, intent(inout), dimension(:,:,:):: kcf
    integer, allocatable, intent(in) :: ptype(:)
    integer, intent(in) :: n
    real, intent(in)    :: sk, gamma
    integer             :: dim, ttp, dtp, ivr

    call getdim(dim)
    call get_tasktype(ttp)
    call get_difftype(dtp)
    call ginitvar(ivr)

    dcf(:,:) = 0.

    select case (ttp)
    case (1)
      ! hydroshock
      call findneighboursKDT(ptype, pos, h)
      call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
      if ( dim > 1 ) then
        call system_clock(start)
        call periodic1v2(den, 20)
        call periodic1v2(h,   20)
        call periodic1v2(om,  20)
        call system_clock(finish)
        call addTime(' bc', finish - start)
      end if
      call eos_adiabatic(n, den, uei, prs, c, gamma)
      if ( dim > 1 ) then
        call system_clock(start)
        call periodic1v2(prs, 20)
        call periodic1v2(c,   20)
        call system_clock(finish)
        call addTime(' bc', finish - start)
      end if
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      if ( dim > 1 ) then
        call system_clock(start)
        call periodic3v2(acc, 20)
        call periodic1v2(due, 20)
        call periodic1v2(dh,  20)
        call system_clock(finish)
        call addTime(' bc', finish - start)
      end if
    case (2)
      ! infslb
      call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
    case (3)
      ! heatconduction
      call findneighboursKDT(ptype, pos, h)

      call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
      call system_clock(start)
      call periodic1v2(den, dim*10)
      call periodic1v2(h,   dim*10)
      call periodic1v2(om,  dim*10)
      call periodic3v2(dcf, dim*10)
      call system_clock(finish)
      call addTime(' BC', finish - start)

      ! symm-diff case for two first derivatives
      ! call c15(pos, mas, h, den, cf, om, dcf)
      ! call system_clock(start)
      ! call periodic3v2(dcf, dim*10)
      ! call system_clock(finish)
      ! call addTime(' BC', finish - start)

      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      if (ivr == 3) then
        if (dim > 1) then
          call system_clock(start)
          call periodic3v2(dcf, dim*10)
          call system_clock(finish)
          call addTime(' BC', finish - start)
        end if
      end if
    case(4)
    case(5)
      ! 'diff-laplace'
      call findneighboursN2plusStatic(ptype, pos, h)
      select case(dtp)
      case(1)
        call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      case(2)
        call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
        call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      case default
        print *, 'Diff type is not set in iterator'
        stop
      end select
    case(6)
      ! 'diff-graddiv'
      call findneighboursN2plusStatic(ptype, pos, h)
      call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
    case(7,8)
    case(10)
      ! diff-artvisc
      ! call findneighboursN2plus(ptype, pos, h)
      call findneighboursKDT(ptype, pos, h)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
    case(9)
      ! soundwave
      call findneighboursKDT(ptype, pos, h)
      call c1(ptype, pos, mas, sk, h, den, om, cf, dcf, kcf)
      call periodic1indims(den, dim)
      call periodic1indims(h, dim)

      call eos_isothermal(den, c(1), prs)

      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)

      call periodic3(acc, 00, dim)
    case default
      print *, 'Task type was not defined in iterator.f90: line 140.'
      stop
    end select
  end subroutine iterate

  subroutine periodic1indims(arr, dim)
    real, allocatable, intent(inout) :: arr(:)
    integer, intent(in) :: dim

    call periodic1(arr, 1)
    if (dim > 1) then
      call periodic1(arr, 2)
      if (dim == 3) then
        call periodic1(arr, 3)
      end if
    end if
  end subroutine
end module
