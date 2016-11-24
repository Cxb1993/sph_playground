module iterator
  use eos
  use circuit1
  use circuit2
  use BC
  use kernel

 implicit none

 public :: iterate

 private

contains
  subroutine iterate(n, sk, gamma, pos, vel, acc, &
                    mas, den, h, dh, om, prs, c, uei, due, cf, dcf, kcf)
    integer, intent(in) :: n
    real, intent(in)    :: sk, gamma
    real, intent(inout) :: pos(3,n), vel(3,n), acc(3,n)
    real, intent(inout) :: mas(n), den(n), h(n), dh(n), prs(n), c(n), uei(n), due(n), om(n)
    real, intent(inout) :: cf(n), dcf(n), kcf(n)
    character (len=40)  :: t
    integer             :: dim

    call get_dim(dim)
    call get_tasktype(t)

    select case (t)
    case ('hydroshock')
      call c1(n, pos, mas, sk, h, den, om)
      call eos_adiabatic(n, den, uei, prs, c, gamma)
      call c2(n, c, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      if (dim.gt.0) then
        call fixed3(acc, 11, 1, 0.)
        call fixed3(acc, 12, 1, 0.)
        if (dim.gt.1) then
          call fixed3(acc, 21, 2, 0.)
          call fixed3(acc, 22, 2, 0.)
        end if
      end if
    case ('infslb')
      call c1(n, pos, mas, sk, h, den, om)
      call c2(n, c, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
    case ('hc-sinx')
      call c1(n, pos, mas, sk, h, den, om)
      ! print *, den
      call periodic1(den, 1)
      call periodic1(h, 1)
      ! read*
      if (dim > 1) then
        call periodic1(den, 2)
        call periodic1(h, 2)
        if (dim == 3) then
          call periodic1(den, 3)
          call periodic1(h, 3)
        end if
      end if
      ! print *, den
      call c2(n, c, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      call periodic1(due, 1)
      if (dim > 1) then
        call periodic1(due, 2)
        if (dim == 3) then
          call periodic1(due, 3)
        end if
      end if
    end select
  end subroutine iterate
end module iterator
