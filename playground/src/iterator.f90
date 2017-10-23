module iterator
  use eos
  use circuit1
  use timing,           only: addTime
  use circuit2,         only: c2
  use BC
  use state,            only: get_difftype,&
                              getdim,&
                              get_tasktype, &
                              ginitvar
  use neighboursearch,  only: findneighboursN2plus, &
                              findneighboursN2, &
                              findneighboursKDT

 implicit none

 public :: iterate

 private
 save
 integer(8)  :: start=0, finish=0

contains
  subroutine iterate(n, sk, gamma, ptype, pos, vel, acc, &
                    mas, den, h, dh, om, prs, c, uei, due, cf, dcf, kcf, dfdx)
    real, allocatable, intent(inout), dimension(:,:)  :: pos, vel, acc, cf, dcf
    real, allocatable, intent(inout), dimension(:)    :: mas, den, dh, prs, c, uei, due, om, h
    real, allocatable, intent(inout), dimension(:,:,:):: dfdx, kcf
    integer, allocatable, intent(in) :: ptype(:)
    integer, intent(in) :: n
    real, intent(in)    :: sk, gamma
    integer             :: dim, ttp, dtp, ivr

    call getdim(dim)
    call get_tasktype(ttp)
    call get_difftype(dtp)
    call ginitvar(ivr)

    select case (ttp)
    case (1)
      ! hydroshock
      call findneighboursKDT(ptype, pos, h)
      call c1(pos, mas, vel, sk, h, den, om, dfdx)
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
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      if ( dim > 1 ) then
        call system_clock(start)
        call periodic3v2(acc, 20)
        call periodic1v2(due, 20)
        call periodic1v2(dh, 20)
        call system_clock(finish)
        call addTime(' bc', finish - start)
      end if
    case (2)
      ! infslb
      call c1(pos, mas, vel, sk, h, den, om, dfdx)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
    case (3)
      ! heatconduction
      call findneighboursN2plus(ptype, pos, h)
      call c1(pos, mas, vel, sk, h, den, om, dfdx)
      ! call periodic1indims(den, dim)
      ! call periodic1indims(h, dim)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      if (ivr == 3) then
        if (dim > 1) then
          call system_clock(start)
          call periodic3(cf, 20, dim)
          if ( dim == 3) then
            call periodic3(cf, 30, dim)
          end if
          call system_clock(finish)
          call addTime(' BC', finish - start)
        end if
      end if
      ! call fixed1(dcf, 11, 0.)
      ! call fixed1(dcf, 12, 0.)
      ! if (dim > 1) then
      !   call fixed1(dcf, 21, 0.)
      !   call fixed1(dcf, 22, 0.)
      !   if (dim == 3) then
      !     call fixed1(dcf, 31, 0.)
      !     call fixed1(dcf, 32, 0.)
      !   end if
      ! end if
    case(4)
    case(5)
      ! 'diff-laplace'
      call findneighboursN2plus(ptype, pos, h)
      select case(dtp)
      case(1)
        call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      case(2)
        call c1(pos, mas, vel, sk, h, den, om, dfdx)
        call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      case default
        print *, 'Diff type is not set in iterator'
        stop
      end select
    case(6)
      ! 'diff-graddiv'
      call findneighboursN2plus(ptype, pos, h)
      call c1(pos, mas, vel, sk, h, den, om, dfdx)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
    case(7,8)
    case(10)
      ! diff-artvisc
      call findneighboursN2plus(ptype, pos, h)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
    case(9)
      ! soundwave
      call findneighboursN2(ptype, pos, h)
      call c1(pos, mas, vel, sk, h, den, om, dfdx)
      call periodic1indims(den, dim)
      call periodic1indims(h, dim)

      call eos_isothermal(den, c(1), prs)

      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)

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
