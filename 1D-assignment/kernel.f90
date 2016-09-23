module kernel
  implicit none

  public :: get_kernel

  private

 contains

  subroutine get_kernel(q, w, dw)
    real, intent(in) :: q
    real             :: w, dw

    if (q >= 2) then
      w  = 0.
      dw = 0.
    else if (q >= 1) then
      w  = 0.25 * ((2. - q) ** 3)
      dw = (- 0.75 * ((2. - q) ** 2)) / q
    else if (q >= 0) then
      w  = 1. - 1.5 * q**2 + 0.75 * q ** 3
      dw = -3. + 2.25 * q
      ! dw = (- 0.75 * ((2. - q) ** 2) + 3.0 * ((1. - q) ** 2)) / q
    else
      print *, 'something went wrong, q =', q
      w  = 0.
      dw = 0.
    end if
    ! w = -0.06*q**2 - 0.2*q + 0.04*(7.07*q - 7.9)*atan(7.07*q - 7.9) - 0.02*log((7.07*q - 7.9)**2 + 1) + 0.36
    ! dw = -0.13*q + 0.32*atan(7.07*q - 7.9) - 0.2
     w = 2./3. * w
    dw = 2./3. * dw
  end subroutine get_kernel

end module kernel
