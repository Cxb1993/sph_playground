module internal
  use eos
  use circuit1
  use circuit2
  use BC
  use kernel

 implicit none

 public :: derivs

 private

contains
  subroutine derivs(n, sk, gamma, pos, vel, acc, &
                    mas, den, h, dh, om, prs, c, uei, due, cf, dcf, kcf)
    integer, intent(in) :: n
    real, intent(in)    :: sk, gamma
    real, intent(inout) :: pos(n,3), vel(n,3), acc(n,3)
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
        call set_fixed3(acc, 11, 1, 0.)
        call set_fixed3(acc, 12, 1, 0.)
        if (dim.gt.1) then
          call set_fixed3(acc, 21, 2, 0.)
          call set_fixed3(acc, 22, 2, 0.)
        end if
      end if
    case ('heatslab')
      call c1(n, pos, mas, sk, h, den, om)
      call c2(n, c, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, cf, dcf, kcf)
      ! if (dim.gt.0) then
      !   call set_fixed1(due, 11, 0.)
      !   call set_fixed1(due, 12, 0.)
      !   if (dim.gt.1) then
      !     call set_fixed1(due, 21, 0.)
      !     call set_fixed1(due, 22, 0.)
      !   end if
      ! end if
    end select
  end subroutine derivs
end module internal
