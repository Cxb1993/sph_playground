module kernel
  use const
  ! use cubic
  use n2movedgaus
  ! use n2ext
  use quintic
  ! use gaus
  ! use sinc
  ! use external

  implicit none

  public :: set_dim, get_nw, get_dw_dh, get_w, get_dim,             &
            set_tasktype, get_tasktype, set_kerntype, get_kerntype, &
            get_n2w, get_krad, get_jacobian, GradDivW, PureKernel!, get_n2y !, get_dphi_dh,

  private
  save
    integer :: dim = 1
    integer :: ttype, ktype

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
     case('diff-laplace')
       ttype = 5
     case('diff-graddiv')
       ttype = 6
     case('chi-laplace')
       ttype = 7
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
    !  call calc_params()
   end subroutine set_kerntype

   subroutine get_kerntype(ott)
     integer, intent(out) :: ott
     ott = ktype
   end subroutine get_kerntype

   subroutine get_kernelname(kname)
     character (len=*), intent(out) :: kname
     kname = kernelname
   end subroutine get_kernelname


   subroutine get_krad(kr)
     real, intent(out) :: kr
     kr = krad
   end subroutine get_krad

  subroutine get_w(r, h, w)
    real, intent(in)  :: r, h
    real, intent(out) :: w
    real              :: f

    call kf(r, h, f)
    w = knorm(dim) * f / h ** dim
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
    real              :: df, ddf, r2, dr

    call kddf(r, h, ddf)
    call kdf(r, h, df)
    ! n2w = knorm(dim)*(ddf + (dim/r - 1) * df)/h**(dim+2)
    n2w = knorm(dim)*(ddf + (dim - 1) * df)/h**(dim+2)
  end subroutine get_on2w

  subroutine get_on2iw(r, h, n2w)
    real, intent(in)    :: r(3), h
    real, intent(out)   :: n2w(3)
    real                :: r2, dr, km(3), df, ddf

    r2 = dot_product(r,r)
    dr = sqrt(r2)
    km(:) = r(:)*r(:)/r2
    call kddf(dr, h, ddf)
    call kdf(dr, h, df)

    n2w(1:dim) = knorm(dim)*(ddf*km(1:dim) + (1 - km(1:dim)) * df)/h**(dim+2)
    ! print *, n2w, ddf*km(1:dim) + (1 - km(1:dim))
    ! read *
  end subroutine get_on2iw

  subroutine get_Fabiw(r, h, Fabw)
    real, intent(in)    :: r(3), h
    real, intent(out)   :: Fabw(3)
    real                :: nw(3)

    call get_nw(r, h, nw)
    Fabw(:) = -2. * r(:) * nw(:)/dot_product(r,r)
  end subroutine get_Fabiw

  subroutine get_n2w(r, h, n2w)
    real, intent(in)  :: r(3), h
    real, intent(out) :: n2w

    if (ktype == 1) then
      call get_on2w(sqrt(dot_product(r,r)), h, n2w)
    else if (ktype == 2) then
      call get_Fab(r, h, n2w)
    end if
  end subroutine get_n2w

  subroutine get_jacobian(r, h, J)
    real, intent(in)  :: r(3), h
    real, intent(out) :: J(3,3)
    real              :: r2, dr, df, ddf

    r2 = dot_product(r,r)
    dr = sqrt(r2)

    call kddf(dr, h, ddf)
    call kdf(dr, h, df)
    J(1,1) = knorm(dim)*(ddf*r(1)*r(1)/r2 + df*(1/dr - r(1)*r(1)/r2))/h**(dim+2) ! d2/dx2   ! Wxx
    J(1,2) = knorm(dim)*(ddf*r(2)*r(1)/r2 - df*r(2)*r(1)/r2)/h**(dim+2)          ! d2/dydx  ! Wxy
    J(1,3) = knorm(dim)*(ddf*r(3)*r(1)/r2 - df*r(3)*r(1)/r2)/h**(dim+2)          ! d2/dzdx  ! Wxz

    J(2,1) = knorm(dim)*(ddf*r(1)*r(2)/r2 - df*r(1)*r(2)/r2)/h**(dim+2)          ! d2/dxdy  ! Wyx
    J(2,2) = knorm(dim)*(ddf*r(2)*r(2)/r2 + df*(1/dr - r(2)*r(2)/r2))/h**(dim+2) ! d2/dy2   ! Wyy
    J(2,3) = knorm(dim)*(ddf*r(3)*r(2)/r2 - df*r(3)*r(2)/r2)/h**(dim+2)          ! d2/dxdz  ! Wyz

    J(3,1) = knorm(dim)*(ddf*r(1)*r(3)/r2 - df*r(1)*r(3)/r2)/h**(dim+2)          ! d2/dxdz  ! Wzx
    J(3,2) = knorm(dim)*(ddf*r(2)*r(3)/r2 - df*r(2)*r(3)/r2)/h**(dim+2)          ! d2/dydz  ! Wzy
    J(3,3) = knorm(dim)*(ddf*r(3)*r(3)/r2 + df*(1/dr - r(3)*r(3)/r2))/h**(dim+2) ! d2/dz2   ! Wzz
  end subroutine get_jacobian

  subroutine GradDivW(r, h, n2w)
    real, intent(in)    :: r(3), h
    real, intent(out)   :: n2w(3)

    if (ktype == 1) then
      call get_on2iw(r, h, n2w)
    else if (ktype == 2) then
      call get_Fabiw(r, h, n2w)
    end if
  end subroutine GradDivW

  subroutine PureKernel(dr ,h, df, ddf)
    real, intent(in)  :: dr, h
    real, intent(out) :: df, ddf

    call kdf(dr, h, df)
    call kddf(dr, h, ddf)
    ddf = knorm(dim)*ddf/h**(dim+2)
    df  = knorm(dim)* df/h**(dim+2)
  end subroutine PureKernel

! ---------!
! Y kernel !--------------------------------------------------------------------
!----------!

  subroutine get_nY(rab, h, ny)
    real, intent(in)  :: rab(3), h
    real, intent(out) :: ny(3)
    real              :: df

    call n2df(sqrt(dot_product(rab(:),rab(:))), h, df)

    ny(:) = n2C(dim,ktype) * df * rab(:) / h**(dim+2)
  end subroutine get_nY

  subroutine get_FabY(r, h, FabY)
    real, intent(in)  :: r(3), h
    real, intent(out) :: FabY
    real              :: nY(3)

    call get_nY(r, h, nY)
    FabY = -2. * dot_product(r,nY)/dot_product(r,r)
  end subroutine get_FabY

  subroutine get_on2Y(r, h, n2Y)
    real, intent(in)  :: r, h
    real, intent(out) :: n2Y
    real              :: df, ddf

    call n2ddf(r, h, ddf)
    call n2df(r, h, df)
    n2Y = n2C(dim,ktype)*(ddf + (dim - 1) * df)/h**(dim+2)
  end subroutine get_on2Y

  subroutine get_on2iY(r, h, n2Y)
    real, intent(in)    :: r(3), h
    real, intent(out)   :: n2Y(3)
    real                :: r2, dr, km(3), df, ddf

    r2 = dot_product(r,r)
    dr = sqrt(r2)
    km(:) = r(:)*r(:)/r2
    call n2ddf(dr, h, ddf)
    call n2df(dr, h, df)
    n2Y(1:dim) = n2C(dim,ktype)*(ddf*km(1:dim) + (1 - km(1:dim)) * df)/h**(dim+2)
  end subroutine get_on2iY

  subroutine get_FabiY(r, h, FabY)
    real, intent(in)    :: r(3), h
    real, intent(out)   :: FabY(3)
    real                :: nY(3)

    call get_nY(r, h, nY)
    FabY(:) = -2. * r(:) * nY(:)/dot_product(r,r)
  end subroutine get_FabiY

  ! subroutine get_n2w(r, h, n2Y)
  !   real, intent(in)  :: r(3), h
  !   real, intent(out) :: n2Y
  !
  !   if (ktype == 1) then
  !     call get_on2Y(sqrt(dot_product(r,r)), h, n2Y)
  !   else if (ktype == 2) then
  !     call get_FabY(r, h, n2Y)
  !   end if
  ! end subroutine get_n2w
  !
  ! subroutine GradDivW(r, h, n2Y)
  !   real, intent(in)    :: r(3), h
  !   real, intent(out)   :: n2Y(3)
  !
  !   if (ktype == 1) then
  !     call get_on2iY(r, h, n2Y)
  !   else if (ktype == 2) then
  !     call get_FabiY(r, h, n2Y)
  !   end if
  ! end subroutine GradDivW

end module kernel