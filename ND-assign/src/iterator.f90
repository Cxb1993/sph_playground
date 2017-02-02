module iterator
  use eos
  use circuit1
  use circuit2,        only:  c2
  use BC
  use kernel
  use neighboursearch, only:  findneighbours

 implicit none

 public :: iterate

 private

contains
  subroutine iterate(n, sk, gamma, ptype, pos, vel, acc, &
                    mas, den, h, dh, om, prs, c, uei, due, cf, dcf, kcf)
    real, allocatable, intent(inout) :: pos(:,:), vel(:,:), acc(:,:), mas(:), den(:), &
                                        dh(:), prs(:), c(:), uei(:), due(:), om(:), &
                                        cf(:), dcf(:), kcf(:),  h(:)
    integer, allocatable, intent(in) :: ptype(:)
    integer, intent(in) :: n
    real, intent(in)    :: sk, gamma
    integer             :: t
    integer             :: dim
    real                :: start, finish

    call get_dim(dim)
    call get_tasktype(t)

    select case (t)
    case (1)
      ! hydroshock
      call c1(n, pos, mas, sk, h, den, om)
      call eos_adiabatic(n, den, uei, prs, c, cf, gamma)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
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
      call c1(n, pos, mas, sk, h, den, om)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
    case (3)
      ! hc-sinx
      call c1(n, pos, mas, sk, h, den, om)
      call periodic1indims(den, dim)
      call periodic1indims(h, dim)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
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
      call c1(n, pos, mas, sk, h, den, om)
      call periodic1indims(den, dim)
      call periodic1indims(h, dim)
      ! print *, '----- 1'
      call eos_adiabatic(n, den, uei, prs, c, cf, gamma)
      call periodic1indims(prs, dim)
      call periodic1indims(c, dim)
      ! print *, '----- 2'
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
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
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      call periodic3(acc, 00, dim)
    case(6)
      ! 'diff-graddiv'
      call cpu_time(start)
      call findneighbours(pos, ptype, h)
      call cpu_time(finish)
      print '("Time in NS = ",f9.4," seconds.")',finish-start
      call cpu_time(start)
      call c2(c, ptype, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      call cpu_time(finish)
      print '("Time in C2 = ",f9.4," seconds.")',finish-start
      call cpu_time(start)
      call periodic3(acc, 00, dim)
      print '("Time in BC = ",f9.4," seconds.")',finish-start
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
