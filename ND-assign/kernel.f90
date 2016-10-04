module kernel
  implicit none

  public :: set_dim, get_nabla_w, get_dw_dh
  real, save :: dim = 1

  private

 contains
   subroutine set_dim(d)
     real, intent(in) :: d
     dim = d
   end subroutine set_dim

  subroutine get_kernel(r, h, w, dw)
    real, intent(in)  :: r, h
    real, intent(out) :: w, dw
    real              :: q

    q = r / h
    if (q >= 2.) then
      w  = 0.
      dw = 0.
    else if (q >= 1.) then
      w  = 0.25 * ((2. - q) ** 3)
      dw = (- 0.75 * ((2. - q) ** 2)) / q
    else if (q >= 0.) then
      w  = 1. - 1.5 * q**2 + 0.75 * q ** 3
      dw = -3. + 2.25 * q
    else
      print *, 'something went wrong, q =', q
      stop
    end if
     w = 2./3. * w
    dw = 2./3. * dw
  end subroutine get_kernel

  subroutine get_nabla_w(rab, h, w, nw)
    real, intent(in)  :: rab(3), h
    real, intent(out) :: w, nw(3)
    real              :: f, df

    call get_kernel(sqrt(dot_product(rab(:),rab(:))), h, f, df)
    w = f / h
    nw(:) = df * rab(:) / h**(dim+2)
  end subroutine get_nabla_w

  subroutine get_dw_dh(r, h, w, dwdh)
    real, intent(in)  :: r, h
    real, intent(out) :: w, dwdh
    real              :: f, df

    call get_kernel(r, h, f, df)
    w = f / h
    dwdh = - (dim * f + r * df / h) / h ** (dim + 1)
  end subroutine get_dw_dh
end module kernel
