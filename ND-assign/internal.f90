module internal
  use eos
  use purehydro_setup
  use purehydro_circuit1
  use purehydro_circuit2
  use tempr_setup
  use tempr_circuit1
  use tempr_circuit2

 implicit none

 public :: derivs

 private

contains
  subroutine derivs(t, n, sos, sk, gamma, pos, vel, acc, mas, den, h, dh, om, prs, c, uei, due, cf, dcf, kcf)
    integer, intent(in)           :: n
    real, intent(in)              :: sos, sk, gamma
    real, intent(inout)           :: pos(n,3), vel(n,3), acc(n,3)
    real, intent(inout)           :: mas(n), den(n), h(n), dh(n), prs(n), c(n), uei(n), due(n), om(n)
    real, intent(inout)           :: cf(n), dcf(n), kcf(n)
    character (len=*), intent(in) :: t

    select case (t)
    case ('periodic')
      call purehydro_c1(n, pos, mas, sk, h, den, om)
      call set_periodic(n, den)
      call eos_isothermal(n, den, prs, sos)
      call set_periodic(n, prs)
      call purehydro_c2(n, c, pos, vel, acc, mas, den, h, om, prs, uei, due, dh)
      call set_periodic(n, acc)
    case ('purehydroshock')
      call purehydro_c1(n, pos, mas, sk, h, den, om)
      call eos_adiabatic(n, den, uei, prs, c, gamma)
      call purehydro_c2(n, c, pos, vel, acc, mas, den, h, om, prs, uei, due, dh)
      call purehydro_set_fixed3(acc)
      call purehydro_set_fixed1(due)
    case ('temperhomog01')
      call tempr_solid_c1(n, pos, mas, sk, h, den)
      call tempr_solid_c1(n, pos, mas, sk, h, den)
      call tempr_solid_c1(n, pos, mas, sk, h, den)
      call tempr_solid_c2(n, pos, mas, den, h, uei, due, kcf)
    end select
  end subroutine derivs
end module internal
