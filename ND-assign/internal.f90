module internal
  use eos
  use circuit1
  use circuit2
  use BC

 implicit none

 public :: derivs

 private

contains
  subroutine derivs(n, sos, sk, gamma, pos, vel, acc, mas, den, h, dh, om, prs, c, uei, due, cf, dcf, kcf)
    integer, intent(in) :: n
    real, intent(in)    :: sos, sk, gamma
    real, intent(inout) :: pos(n,3), vel(n,3), acc(n,3)
    real, intent(inout) :: mas(n), den(n), h(n), dh(n), prs(n), c(n), uei(n), due(n), om(n)
    real, intent(inout) :: cf(n), dcf(n), kcf(n)
    character (len=40)   :: t

    call get_tasktype(t)
    select case (t)
    case ('hydroshock')
      call make_c1(n, pos, mas, sk, h, den, om)
      call eos_adiabatic(n, den, uei, prs, c, gamma)
      call make_c2(n, c, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, kcf)
      call set_fixed3(acc)
      call set_fixed1(due)
    case ('temperhomog01')
      call make_c1(n, pos, mas, sk, h, den, om)
      call make_c2(n, c, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, kcf)
    end select
  end subroutine derivs
end module internal
