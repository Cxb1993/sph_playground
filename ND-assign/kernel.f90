module kernel
  implicit none

  public :: set_dim, get_nabla_w, get_dw_dh, get_w, get_dim
  integer, save   :: dim = 1
  real, parameter :: pi = 4.*atan(1.)

  private

 contains
   subroutine set_dim(d)
     integer, intent(in) :: d
     dim = d
   end subroutine set_dim

   subroutine get_dim(d)
     integer, intent(out) :: d
     d = dim
   end subroutine get_dim

  subroutine get_kernel_f(r, h, f)
    real, intent(in)  :: r, h
    real, intent(out) :: f
    real              :: q

    q = r / h
    if (q >= 2.) then
      f  = 0.
    else if (q >= 1.) then
      f  = 0.25 * ((2. - q) ** 3)
    else if (q >= 0.) then
      f  = 1. - 1.5 * q**2 + 0.75 * q ** 3
    else
      print *, 'something went wrong, q =', q
      stop
    end if
    select case(dim)
      case(1)
        f = 2./3. * f
      case(2)
        f = 10./(7. * pi) * f
      case(3)
        f = f / pi
    end select
  end subroutine get_kernel_f

  subroutine get_kernel_df(r, h, df)
    real, intent(in)  :: r, h
    real, intent(out) :: df
    real              :: q

    q = r / h
    if (q >= 2.) then
      df = 0.
    else if (q >= 1.) then
      df = (- 0.75 * ((2. - q) ** 2)) / q
    else if (q >= 0.) then
      df = -3. + 2.25 * q
    else
      print *, 'something went wrong, q =', q
      stop
    end if
    select case(dim)
      case(1)
        df = 2./3. * df
      case(2)
        df = 10./(7. * pi) * df
      case(3)
        df = 1./pi * df
    end select
  end subroutine get_kernel_df

  subroutine get_nabla_w(rab, h, nw)
    real, intent(in)  :: rab(3), h
    real, intent(out) :: nw(3)
    real              :: df
    integer :: i

    call get_kernel_df(sqrt(dot_product(rab(:),rab(:))), h, df)

    nw(:) = df * rab(:) / h**(dim+2)
  end subroutine get_nabla_w

  subroutine get_dw_dh(r, h, dwdh)
    real, intent(in)  :: r, h
    real, intent(out) :: dwdh
    real              :: f, df

    call get_kernel_f(r, h, f)
    call get_kernel_df(r, h, df)

    dwdh = - (dim * f + r * df / h) / h ** (dim + 1)
  end subroutine get_dw_dh

  subroutine get_w(r, h, w)
    real, intent(in)  :: r, h
    real, intent(out) :: w
    real              :: f

    call get_kernel_f(r, h, f)

    w = f / h ** dim
  end subroutine get_w
end module kernel
