module kernel
  implicit none

  public :: set_dim, get_nabla_w, get_dw_dh, get_w, get_dim, get_dphi_dh, get_n2y
  integer, save   :: dim = 1
  real, parameter :: pi = 4.*atan(1.)

  private

 contains
   subroutine set_dim(d)
     integer, intent(in) :: d
     dim = d
    !  print *, 'Dim in set call: ', dim
    !  read *
   end subroutine set_dim

   subroutine get_dim(d)
     integer, intent(out) :: d
     d = dim
    !  print *, 'Dim in get call: ', d
    !  read *
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

  subroutine get_kernel_phi(r, h, phi)
    real, intent(in)  :: r, h
    real, intent(out) :: phi

    phi = 1./h * (1. + (r/h)**2)**(1./2.)
  end subroutine get_kernel_phi

  subroutine get_w(r, h, w)
    real, intent(in)  :: r, h
    real, intent(out) :: w
    real              :: f

    call get_kernel_f(r, h, f)

    w = f / h ** dim
  end subroutine get_w

  subroutine get_nabla_w(rab, h, nw)
    real, intent(in)  :: rab(3), h
    real, intent(out) :: nw(3)
    real              :: df

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

  subroutine get_dphi_dh(r, h, dphidh)
    real, intent(in)  :: r, h
    real, intent(out) :: dphidh
    real              :: phi

    call get_kernel_phi(r, h, phi)
    dphidh = - 1./h * phi - r**2/(h**5 * phi)
  end subroutine get_dphi_dh

  subroutine get_n2y(rab, h, n2y)
    real, intent(in)  :: rab(3), h
    real, intent(out) :: n2y
    real              :: df, Fab, r

    r = sqrt(dot_product(rab, rab))

    call get_kernel_df(r, h, df)
    Fab = df / h**(dim+2)
    n2y = - 2 * Fab  / r
  end subroutine get_n2y
end module kernel
