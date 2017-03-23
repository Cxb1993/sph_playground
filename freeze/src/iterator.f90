module iterator
  use eos
  use circuit1
  use circuit2,         only:  c2
  use BC
  use kernel,           only: get_difftype,&
                              get_dim,&
                              get_tasktype
  use neighboursearch,  only: findneighbours

 implicit none

 public :: iterate

 private

contains
  subroutine iterate(n, sk, gamma, ptype, pos, vel, acc, &
                    mas, den, h, dh, om, prs, c, uei, due, cf, dcf, kcf, dfdx)
    real, allocatable, intent(inout) :: pos(:,:), vel(:,:), acc(:,:), mas(:), den(:), &
                                        dh(:), prs(:), c(:), uei(:), due(:), om(:), &
                                        cf(:), dcf(:), kcf(:),  h(:), dfdx(:,:,:)
    integer, allocatable, intent(in) :: ptype(:)
    integer, intent(in) :: n
    real, intent(in)    :: sk, gamma
    integer             :: dim, ttp, dtp

    call get_dim(dim)
    call get_tasktype(ttp)
    call get_difftype(dtp)

    select case (ttp)
    case (1)
      ! hydroshock
      call findneighbours(ptype, pos, h)
      call c1(ptype, pos, mas, vel, sk, h, den, om, dfdx)
      call eos_adiabatic(n, den, uei, prs, c, cf, gamma)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      if (dim.gt.0) then
        call fixed3(acc, 11, 1, 0.)
        call fixed3(acc, 12, 1, 0.)
        if (dim.gt.1) then
          call fixed3(acc, 21, 2, 0.)
          call fixed3(acc, 22, 2, 0.)
        end if
      end if
    case (2)
      ! infslb
      call findneighbours(ptype, pos, h)
      call c1(ptype, pos, mas, vel, sk, h, den, om, dfdx)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
    case (3)
      ! hc-sinx
      call findneighbours(ptype, pos, h)
      call c1(ptype, pos, mas, vel, sk, h, den, om, dfdx)
      call periodic1indims(den, dim)
      call periodic1indims(h, dim)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      call periodic1indims(due, dim)
      ! call periodic1(due, 1)
      ! if (dim > 1) then
      !   call periodic1(due, 2)
      !   if (dim == 3) then
      !     call periodic1(due, 3)
      !   end if
      ! end if
    case(4)
      ! print *, pos
      ! print *, den
      ! print *, h
      ! read *
      ! print *, '----- 0'
      call findneighbours(ptype, pos, h)
      call c1(ptype, pos, mas, vel, sk, h, den, om, dfdx)
      call periodic1indims(den, dim)
      call periodic1indims(h, dim)
      ! print *, '----- 1'
      call eos_adiabatic(n, den, uei, prs, c, cf, gamma)
      call periodic1indims(prs, dim)
      call periodic1indims(c, dim)
      ! print *, '----- 2'
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      ! call periodic3(acc, 10, 0)
      ! call periodic3(due, 10, 0)
      ! call periodic3(dcf, 10, 0)
      ! print *, '----- 3'
      ! if (dim.gt.0) then
      !   call fixed3(acc, 11, 1, 0.)
      !   call fixed3(acc, 12, 1, 0.)
      !   if (dim.gt.1) then
      !     call fixed3(acc, 21, 2, 0.)
      !     call fixed3(acc, 22, 2, 0.)
      !   end if
      ! end if
    case(5)
      ! 'diff-laplace'
      call findneighbours(ptype, pos, h)
      select case(dtp)
      case(1)
        call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      case(2)
        call c1(ptype, pos, mas, vel, sk, h, den, om, dfdx)
        call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      case default
        print *, 'Diff type is not set in iterator'
        stop
      end select
    case(6)
      ! 'diff-graddiv'
      call findneighbours(ptype, pos, h)
      select case(dtp)
      case(1)
        call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      case(2)
        call c1(ptype, pos, mas, vel, sk, h, den, om, dfdx)
        ! print*, '======================'
        call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf, dfdx)
      case default
        print *, 'Diff type is not set in iterator'
        stop
      end select
    case(7,8)
    case default
      print *, 'Task type was not defined in iterator'
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
  end subroutine periodic1indims
end module iterator
