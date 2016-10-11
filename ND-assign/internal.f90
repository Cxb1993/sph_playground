module internal
  use eos
  use setup
  use circuit1_mod
  use circuit2_mod

 implicit none

 public :: derivs

 private

contains

  subroutine derivs(t, n, sos, sk, gamma, pos, vel, acc, mas, den, h, dh, om, prs, c, uei, due, f, eps)
    integer, intent(in)           :: n
    real, intent(in)              :: sos, sk, gamma
    real, intent(inout)           :: pos(n,3), vel(n,3), acc(n,3)
    real, intent(inout)           :: mas(n), den(n), h(n), dh(n), prs(n), c(n), uei(n), due(n), om(n), f(n), eps(n)
    character (len=*), intent(in) :: t

    call circuit1(n, pos, mas, sk, h, den, om)
    if (t.eq.'periodic') then
      call set_periodic(n, den)
    end if

    select case (t)
      case ('periodic')
        call eos_isothermal(n, den, prs, sos)
      case ('shock_fixed')
        call eos_adiabatic(n, den, uei, prs, c, gamma)
    end select

    if (t.eq.'periodic') then
      call set_periodic(n, prs)
    end if

    call circuit2(n, c, pos, vel, acc, mas, den, h, om, prs, uei, due, dh, f, eps)

    select case (t)
      case ('periodic')
        call set_periodic(n, acc)
      case ('shock_fixed')
        call set_fixed3(acc)
        call set_fixed1(due)
      end select

  end subroutine derivs

end module internal
