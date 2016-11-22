module kernel
  use cubic
  ! use quintic
  implicit none

  public :: set_dim, get_nw, get_dw_dh, get_w, get_dim, &
            set_tasktype, get_tasktype, set_kerntype, get_kerntype, &
            get_n2w, get_krad !get_n2y, get_dphi_dh,

  private
    integer, save            :: dim = 1
    real, parameter          :: pi = 4.*atan(1.)
    character (len=40), save :: ttype, ktype

 contains
  !  GetterSetter accecc methods
   subroutine set_dim(d)
     integer, intent(in) :: d
     dim = d
   end subroutine set_dim

   subroutine get_dim(d)
     integer, intent(out) :: d
     d = dim
   end subroutine get_dim

   subroutine set_tasktype(itt)
     character (len=*), intent(in) :: itt
     ttype = itt
   end subroutine set_tasktype

   subroutine get_tasktype(ott)
     character (len=*), intent(out) :: ott
     ott = ttype
   end subroutine get_tasktype

   subroutine set_kerntype(itt)
     character (len=*), intent(in) :: itt
     ktype = itt
   end subroutine set_kerntype

   subroutine get_kerntype(ott)
     character (len=*), intent(out) :: ott
     ott = ktype
   end subroutine get_kerntype

   subroutine get_krad(kr)
     real, intent(out) :: kr
     kr = krad
   end subroutine get_krad

  subroutine get_w(r, h, w)
    real, intent(in)  :: r, h
    real, intent(out) :: w
    real              :: f

    call kf(r, h, f)
    f = knorm(dim)* f
    w = f / h ** dim
  end subroutine get_w

  subroutine get_nw(rab, h, nw)
    real, intent(in)  :: rab(3), h
    real, intent(out) :: nw(3)
    real              :: df

    call kdf(sqrt(dot_product(rab(:),rab(:))), h, df)

    nw(:) = knorm(dim) * df * rab(:) / h**(dim+2)
  end subroutine get_nw

  subroutine get_dw_dh(r, h, dwdh)
    real, intent(in)  :: r, h
    real, intent(out) :: dwdh
    real              :: f, df

    call kf(r, h, f)
    call kdf(r, h, df)
    dwdh = - knorm(dim) * (dim * f + r * df / h) / h ** (dim + 1)
  end subroutine get_dw_dh

  subroutine get_Fab(r, h, Fab)
    real, intent(in)  :: r(3), h
    real, intent(out) :: Fab
    real              :: nw(3)

    call get_nw(r, h, nw)
    Fab = -2. * dot_product(r,nw)/dot_product(r,r)
  end subroutine get_Fab

  subroutine get_on2w(r, h, n2w)
    real, intent(in)  :: r, h
    real, intent(out) :: n2w
    real              :: df, ddf

    call kddf(r, h, ddf)
    call kdf(r, h, df)
    n2w = knorm(dim)*(ddf + (dim - 1) * df)/h**(dim+2)
  end subroutine get_on2w

subroutine get_n2w(r, h, n2w)
  real, intent(in)  :: r(3), h
  real, intent(out) :: n2w

  if (ktype == 'n2w') then
    call get_on2w(sqrt(dot_product(r,r)), h, n2w)
  else if (ktype == 'fab') then
    call get_Fab(r, h, n2w)
  else
    print *, 'Kernel: type not chosen'
    stop
  end if
end subroutine get_n2w

  ! subroutine get_kernel_phi(r, h, phi)
  !   real, intent(in)  :: r, h
  !   real, intent(out) :: phi
  !
  !   phi = 1./h * (1. + (r/h)**2)**(1./2.)
  ! end subroutine get_kernel_phi
  !
  ! subroutine get_dphi_dh(r, h, dphidh)
  !   real, intent(in)  :: r, h
  !   real, intent(out) :: dphidh
  !   real              :: phi
  !
  !   call get_kernel_phi(r, h, phi)
  !   dphidh = - 1./h * phi - r**2/(h**5 * phi)
  ! end subroutine get_dphi_dh
end module kernel
