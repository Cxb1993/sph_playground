module kernel
  use const
  use state
  ! use cubic
  ! use n2movedgaus
  ! use n2ext
  ! use n2fromfabcubic
  ! use n2fromWcubic
  ! use quintic
  ! use gaus
  use sinc
  ! use external
  implicit none

  public :: get_nw, get_dw_dh, get_w, setdimkernel, &
            get_n2w, get_krad, get_hessian!, PureKernel!, GradDivW!, get_n2y !, get_dphi_dh,
  save
    integer :: dim
  private
 contains

  pure subroutine get_kernelname(kname)
    character (len=*), intent(out) :: kname
    kname = kernelname
  end subroutine get_kernelname

  pure subroutine get_krad(kr)
    real, intent(out) :: kr
    kr = krad
  end subroutine get_krad

  subroutine setdimkernel(indim)
    integer, intent(in) :: indim
    dim = indim
    call setdimbase(dim)
    call setdim(dim)
  end subroutine
  !
  ! ---------!
  ! W kernel !------------------------------------------------------------------
  !----------!

  pure subroutine get_w(r, h, w)
    real, intent(in)  :: r, h
    real, intent(out) :: w
    real              :: f

    call kf(r, h, f)
    w = wCv * f / h ** dim
  end subroutine get_w

  pure subroutine get_nw(rab, h, nw)
    real, intent(in)  :: rab(3), h
    real, intent(out) :: nw(3)
    real              :: df

    call kdf(sqrt(dot_product(rab(:),rab(:))), h, df)

    nw(:) = wCv * df * rab(:) / h**(dim+2)
  end subroutine get_nw

  pure subroutine get_dw_dh(r, h, dwdh)
    real, intent(in)  :: r, h
    real, intent(out) :: dwdh
    real              :: f, df

    call kf(r, h, f)
    call kdf(r, h, df)
    dwdh = - wCv * (dim * f + r * df / h) / h ** (dim + 1)
  end subroutine get_dw_dh

  pure subroutine get_Fab(r, h, Fab)
    real, intent(in)  :: r(3), h
    real, intent(out) :: Fab
    real              :: nw(3)

    call get_nw(r, h, nw)
    Fab = -2. * dot_product(r,nw)/dot_product(r,r)
  end subroutine get_Fab

  pure subroutine get_on2w(r, h, n2w)
    real, intent(in)  :: r, h
    real, intent(out) :: n2w
    real              :: df, ddf

    call kddf(r, h, ddf)
    call kdf(r, h, df)
    n2w = wCv*(ddf + (dim - 1) * df)/h**(dim+2)
  end subroutine get_on2w

  pure subroutine get_n2w(r, h, n2w)
    real, intent(in)  :: r(3), h
    real, intent(out) :: n2w
    integer :: ktype

    call get_kerntype(ktype)

    if (ktype == 1) then
      call get_on2w(sqrt(dot_product(r,r)), h, n2w)
    else if (ktype == 2) then
      call get_Fab(r, h, n2w)
    end if
  end subroutine get_n2w

  pure subroutine get_hessian(r, h, Hes)
    real, intent(in)  :: r(3), h
    real, intent(out) :: Hes(3,3)
    real              :: r2, dr, df, ddf, fab
    integer :: ktype

    call get_kerntype(ktype)

    if (ktype == 1) then
      r2 = dot_product(r,r)
      dr = sqrt(r2)

      call kddf(dr, h, ddf)
      call kdf(dr, h, df)
      Hes(1,1) = wCv*(ddf*r(1)*r(1)/r2 + df*(1 - r(1)*r(1)/r2))/h**(dim+2)    ! d2/dx2   ! Wxx
      Hes(1,2) = wCv*(ddf*r(2)*r(1)/r2 - df*r(2)*r(1)/r2)/h**(dim+2)          ! d2/dydx  ! Wxy
      Hes(1,3) = wCv*(ddf*r(3)*r(1)/r2 - df*r(3)*r(1)/r2)/h**(dim+2)          ! d2/dzdx  ! Wxz

      Hes(2,1) = wCv*(ddf*r(1)*r(2)/r2 - df*r(1)*r(2)/r2)/h**(dim+2)          ! d2/dxdy  ! Wyx
      Hes(2,2) = wCv*(ddf*r(2)*r(2)/r2 + df*(1 - r(2)*r(2)/r2))/h**(dim+2)    ! d2/dy2   ! Wyy
      Hes(2,3) = wCv*(ddf*r(3)*r(2)/r2 - df*r(3)*r(2)/r2)/h**(dim+2)          ! d2/dxdz  ! Wyz

      Hes(3,1) = wCv*(ddf*r(1)*r(3)/r2 - df*r(1)*r(3)/r2)/h**(dim+2)          ! d2/dxdz  ! Wzx
      Hes(3,2) = wCv*(ddf*r(2)*r(3)/r2 - df*r(2)*r(3)/r2)/h**(dim+2)          ! d2/dydz  ! Wzy
      Hes(3,3) = wCv*(ddf*r(3)*r(3)/r2 + df*(1 - r(3)*r(3)/r2))/h**(dim+2)    ! d2/dz2   ! Wzz
      if ( dim == 1 ) then
        Hes(1,2:3) = 0.
        Hes(2,:) = 0.
        Hes(3,:) = 0.
      elseif ( dim == 2 ) then
        Hes(3,:) = 0.
        Hes(:,3) = 0.
      end if
    elseif ( ktype == 2 ) then
      r2 = dot_product(r,r)
      dr = sqrt(r2)

      call get_Fab(r, h, fab)
      Hes(1,1) = ((dim+2)*r(1)*r(1)/r2-1)*0.5*fab
      Hes(1,2) = (dim+2)*r(1)*r(2)/r2*0.5*fab
      Hes(1,3) = (dim+2)*r(1)*r(3)/r2*0.5*fab

      Hes(2,1) = (dim+2)*r(2)*r(1)/r2*0.5*fab
      Hes(2,2) = ((dim+2)*r(2)*r(2)/r2-1)*0.5*fab
      Hes(2,3) = (dim+2)*r(2)*r(3)/r2*0.5*fab

      Hes(3,1) = (dim+2)*r(3)*r(1)/r2*0.5*fab
      Hes(3,2) = (dim+2)*r(3)*r(2)/r2*0.5*fab
      Hes(3,3) = ((dim+2)*r(3)*r(3)/r2-1)*0.5*fab
      ! H = 2./3. * H
      if ( dim == 1 ) then
        Hes(1,2:3) = 0.
        Hes(2,:) = 0.
        Hes(3,:) = 0.
      elseif ( dim == 2 ) then
        Hes(3,:) = 0.
        Hes(:,3) = 0.
      end if
    end if
  end subroutine get_hessian
end module kernel
