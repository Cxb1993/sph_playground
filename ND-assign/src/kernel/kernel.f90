module kernel
  ! use cubic
  ! use quintic
  use gaus
  ! use external
  implicit none

  public :: set_dim, get_nw, get_dw_dh, get_w, get_dim, &
            set_tasktype, get_tasktype, set_kerntype, get_kerntype, &
            get_n2w, get_n2iw, get_krad !get_n2y, get_dphi_dh,

  private
    integer, save   :: dim = 1
    real, parameter :: pi = 4.*atan(1.)
    integer, save   :: ttype, ktype

 contains
   !
   !-- GetterSetter access
   !
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
     select case(itt)
     case('hydroshock')
       ttype = 1
     case('infslb')
       ttype = 2
     case('hc-sinx')
       ttype = 3
     case('pheva')
       ttype = 4
     case default
       print *, 'Task type not set: ', itt
       stop
     end select
   end subroutine set_tasktype

   subroutine get_tasktype(ott)
     integer, intent(out) :: ott
     ott = ttype
   end subroutine get_tasktype

   subroutine set_kerntype(itt)
     character (len=*), intent(in) :: itt
     select case(itt)
     case('n2w')
       ktype = 1
     case('fab')
       ktype = 2
     case default
       print *, 'Kernel type not set: ', itt
       stop
     end select
     call calc_params()
   end subroutine set_kerntype

   subroutine get_kerntype(ott)
     integer, intent(out) :: ott
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
    f = knorm(dim) * f
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

  subroutine get_on2iw(r, h, n2w, i)
    real, intent(in)    :: r(3), h
    integer, intent(in) :: i
    real, intent(out)   :: n2w
    real                :: r2, dr, km, df, ddf

    r2 = dot_product(r,r)
    dr = sqrt(r2)
    km = r(i)*r(i)/r2
    call kddf(dr, h, ddf)
    call kdf(dr, h, df)
    n2w = knorm(dim)*(ddf*km + (1 - km) * df)/h**(dim+2)
  end subroutine get_on2iw

  subroutine get_Fabi(r, h, Fab, i)
    real, intent(in)    :: r(3), h
    integer, intent(in) :: i
    real, intent(out)   :: Fab
    real                :: nw(3)

    call get_nw(r, h, nw)
    Fab = nw(i)
    nw(:) = 0.
    nw(i) = Fab
    Fab = -2. * dot_product(r,nw)/dot_product(r,r)
  end subroutine get_Fabi

  subroutine get_n2w(r, h, n2w)
    real, intent(in)  :: r(3), h
    real, intent(out) :: n2w

    if (ktype == 1) then
      call get_on2w(sqrt(dot_product(r,r)), h, n2w)
    else if (ktype == 2) then
      call get_Fab(r, h, n2w)
    end if
  end subroutine get_n2w

  subroutine get_n2iw(r, h, n2w, i)
    real, intent(in)    :: r(3), h
    integer, intent(in) :: i
    real, intent(out)   :: n2w

    if (ktype == 1) then
      call get_on2iw(r, h, n2w, i)
    else if (ktype == 2) then
      call get_Fabi(r, h, n2w, i)
    end if
  end subroutine get_n2iw
end module kernel
