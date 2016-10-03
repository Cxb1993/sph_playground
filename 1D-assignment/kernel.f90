module kernel
  implicit none

  public :: get_kernel, get_kernel_dh

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
    else
      print *, 'something went wrong, q =', q
      stop
    end if
    w = 2./3. * w
    dw = 2./3. * dw
  end subroutine get_kernel

  subroutine get_kernel_dh(q, w, dwdh)
    real, intent(in) :: q
    real             :: w, dwdh

    if (q >= 2) then
      w  = 0.
      dwdh = 0.
    else if (q >= 1) then
      w  = 0.25 * ((2. - q) ** 3)
      dwdh = 0.75 * ((2. - q) ** 2) * q
    else if (q >= 0) then
      w  = 1. - 1.5 * q**2 + 0.75 * q ** 3
      dwdh = -3. * q**2 - 2.25 * q ** 3
    else
      print *, 'something went wrong, q =', q
      stop
    end if
    w = 2./3. * w
    dwdh = 2./3. * dwdh
  end subroutine get_kernel_dh
end module kernel
