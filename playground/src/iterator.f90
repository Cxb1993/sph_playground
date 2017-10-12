module iterator
  use eos
  use circuit1
  use timing,           only: addTime
  use circuit2,         only:  c2
  use BC
  use state,            only: get_difftype,&
                              getdim,&
                              get_tasktype, &
                              ginitvar
  use neighboursearch,  only: findneighboursN2plus, &
                              findneighboursN2

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

    ! call findneighboursN2(ptype, pos, h)

    call findneighboursN2plus(ptype, pos, h)

    select case (ttp)
    case (1)
      ! hydroshock
      call c1(pos, mas, vel, sk, h, den, om, dfdx)
      ! call c1a(pos, mas, sk, h, den)
      call eos_adiabatic(n, den, uei, prs, c, gamma)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      if ( dim > 0 ) then
        call fixed3(acc, 11, 1, 0.)
        call fixed3(acc, 12, 1, 0.)
        if ( dim > 1 ) then
          call fixed3(acc, 21, 2, 0.)
          call fixed3(acc, 22, 2, 0.)
        end if
      end if
    case (2)
      ! infslb
      call c1(pos, mas, vel, sk, h, den, om, dfdx)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
    case (3)
      ! heatconduction
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
      select case(dtp)
      case(1)
        call c1(pos, mas, vel, sk, h, den, om, dfdx)
        call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      case(2)
        call c1(pos, mas, vel, sk, h, den, om, dfdx)
        call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      case default
        print *, 'Diff type is not set in iterator'
        stop
      end select
    case(7,8)
    case(9)
      ! soundwave
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
