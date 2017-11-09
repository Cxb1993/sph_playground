module kernel
  use const
  use state
  use base_kernel
  implicit none

  public :: get_nw, get_dw_dh, get_w, setdimkernel, &
            get_nw_cyl, &
            get_n2w, get_krad, get_hessian, getkernelname, get_hessian_rr, getneibnumber
            !, PureKernel!, GradDivW!, get_n2y !, get_dphi_dh,
  save
    integer :: dim
  private
 contains

  pure subroutine getkernelname(kname)
    character (len=*), intent(out) :: kname
    kname = kernelname
  end subroutine

  pure subroutine get_krad(kr)
    real, intent(out) :: kr
    kr = krad
  end subroutine get_krad

  subroutine getneibnumber(nn)
    integer, intent(out) :: nn

    nn = returnneibnum
  end subroutine

  subroutine setdimkernel(indim)
    integer, intent(in) :: indim
    dim = indim
    call setdimbase(dim)
    call setdim(dim)
  end subroutine

  ! ---------!
  ! W kernel !------------------------------------------------------------------
  !----------!

  pure subroutine get_w(r, h, w)
  ! subroutine get_w(r, h, w)
    real, intent(in)  :: r, h
    real, intent(out) :: w
    real              :: f, q

    q = r / h

    call kf(q, f)

    w = wCv * f / h ** dim
  end subroutine

  pure subroutine get_nw(rab, h, nw)
    real, intent(in)  :: rab(3), h
    real, intent(out) :: nw(3)
    real              :: df, q

    q = sqrt(dot_product(rab(:),rab(:))) / h
    call kdf(q, df)

    nw(:) = wCv * df * rab(:) / h**(dim+2) / q
  end subroutine

  pure subroutine get_nw_cyl(ra, rb, h, nw)
    ! get cylindrical nabla kernel for cartesian input
    real, intent(in)  :: ra(3), rb(3), h
    real, intent(out) :: nw(3)
    real              :: df, q, rab(3), cca(3), ccb(3)

    cca(1) = sqrt(ra(1)*ra(1) + ra(2)*ra(2))
    cca(2) = atan(ra(2),ra(1))
    cca(3) = ra(3)

    ccb(1) = sqrt(rb(1)*rb(1) + rb(2)*rb(2))
    ccb(2) = atan(rb(2),rb(1))
    ccb(3) = rb(3)

    rab(:) = ra(:) - rb(:)
    q = sqrt(dot_product(rab(:),rab(:))) / h
    call kdf(q, df)

    nw(1) = wCv / h**(dim+2) / q * df * (cca(1)-ccb(1)+ccb(1)*(1-cos(cca(2)-ccb(2))))
    nw(2) = wCv / h**(dim+2) / q * df * (cca(1)*ccb(1)*sin(cca(2)-ccb(2)))
    nw(3) = wCv / h**(dim+2) / q * df * (cca(3)-ccb(3))
  end subroutine get_nw_cyl

  pure subroutine get_dw_dh(r, h, dwdh)
    real, intent(in)  :: r, h
    real, intent(out) :: dwdh
    real              :: f, df, q
    q = r / h

    call kf(q, f)
    call kdf(q, df)

    dwdh = - wCv / h**(dim + 1) * (dim * f + q * df)
  end subroutine get_dw_dh

  pure subroutine get_FW(r, h, fw)
    real, intent(in)  :: r(3), h
    real, intent(out) :: fw
    real              :: f, dr, q

    dr = sqrt(dot_product(r,r))
    q = dr / h
    call kf(q, f)
    fw = fwc * wCv * f / h ** (dim + 2)
  end subroutine

  pure subroutine get_Fab(r, h, Fab)
    real, intent(in)  :: r(3), h
    real, intent(out) :: Fab
    real              :: q, df

    q = sqrt(dot_product(r,r)) / h
    call kdf(q, df)
    Fab = -2. * wCv * df / h**(dim+2) / q
  end subroutine

  pure subroutine get_on2w(r, h, n2w)
    real, intent(in)  :: r(:), h
    real, intent(out) :: n2w
    real              :: df, ddf, q

    q = sqrt(dot_product(r,r)) / h

    call kddf(q, ddf)
    call kdf(q, df)
    n2w = wCv*(ddf + (dim - 1) * df / q)/h**(dim + 2)
  end subroutine

  ! pure
  pure subroutine get_n2w(r, h, n2w)
    real, intent(in)  :: r(3), h
    real, intent(out) :: n2w
    integer :: ktype

    call getkerntype(ktype)

    if (ktype == 1) then
      call get_on2w(r, h, n2w)
    else if (ktype == 2) then
      call get_Fab(r, h, n2w)
    else if (ktype == 4) then
      call get_FW(r, h, n2w)
    end if
  end subroutine

  pure subroutine get_hessian(r, h, Hes)
  ! subroutine get_hessian(r, h, Hes)
    real, intent(in)  :: r(3), h
    real, intent(out) :: Hes(3,3)
    real              :: r2, dr, df, ddf, fab, q
    real              :: r11, r12, r13, r22, r23, r33, cstart, dfq
    integer :: ktype

    call getkerntype(ktype)

    r2 = dot_product(r,r)

    r11 = r(1)*r(1)/r2
    r12 = r(1)*r(2)/r2
    r13 = r(1)*r(3)/r2
    r22 = r(2)*r(2)/r2
    r23 = r(2)*r(3)/r2
    r33 = r(3)*r(3)/r2
    cstart = wCv/h**(dim+2)
    Hes(:,:) = 0.

    if (ktype == 1) then
      dr = sqrt(r2)

      q = dr / h

      call kddf(q, ddf)
      call kdf(q, df)

      dfq = df/q
      ! Hes(1,1) = wCv*(ddf*r(1)*r(1)/r2 + df*(1 - r(1)*r(1)/r2)/q)/h**(dim+2)    ! d2/dx2   ! Wxx
      ! Hes(1,2) = wCv*(ddf*r(2)*r(1)/r2 - df*r(2)*r(1)/r2/q)/h**(dim+2)          ! d2/dydx  ! Wxy
      ! Hes(1,3) = wCv*(ddf*r(3)*r(1)/r2 - df*r(3)*r(1)/r2/q)/h**(dim+2)          ! d2/dzdx  ! Wxz
      !
      ! Hes(2,1) = wCv*(ddf*r(1)*r(2)/r2 - df*r(1)*r(2)/r2/q)/h**(dim+2)          ! d2/dxdy  ! Wyx
      ! Hes(2,2) = wCv*(ddf*r(2)*r(2)/r2 + df*(1 - r(2)*r(2)/r2)/q)/h**(dim+2)    ! d2/dy2   ! Wyy
      ! Hes(2,3) = wCv*(ddf*r(3)*r(2)/r2 - df*r(3)*r(2)/r2/q)/h**(dim+2)          ! d2/dxdz  ! Wyz
      !
      ! Hes(3,1) = wCv*(ddf*r(1)*r(3)/r2 - df*r(1)*r(3)/r2/q)/h**(dim+2)          ! d2/dxdz  ! Wzx
      ! Hes(3,2) = wCv*(ddf*r(2)*r(3)/r2 - df*r(2)*r(3)/r2/q)/h**(dim+2)          ! d2/dydz  ! Wzy
      ! Hes(3,3) = wCv*(ddf*r(3)*r(3)/r2 + df*(1 - r(3)*r(3)/r2)/q)/h**(dim+2)    ! d2/dz2   ! Wzz

      Hes(1,1) = cstart*(ddf*r11 + dfq*(1 - r11))    ! d2/dx2   ! Wxx

      if ( dim /= 1 ) then
        Hes(1,2) = cstart*(ddf*r12 - dfq*r12)          ! d2/dydx  ! Wxy
        Hes(2,1) = Hes(1,2)                            ! d2/dxdy  ! Wyx
        Hes(2,2) = cstart*(ddf*r22 + dfq*(1 - r22))    ! d2/dy2   ! Wyy
        if ( dim == 3 ) then
          Hes(1,3) = cstart*(ddf*r13 - dfq*r13)          ! d2/dzdx  ! Wxz
          Hes(3,1) = Hes(1,3)                            ! d2/dxdz  ! Wzx
          Hes(2,3) = cstart*(ddf*r23 - dfq*r23)          ! d2/dxdz  ! Wyz
          Hes(3,2) = Hes(2,3)                            ! d2/dydz  ! Wzy
          Hes(3,3) = cstart*(ddf*r33 + dfq*(1 - r33))    ! d2/dz2   ! Wzz
        end if
      end if
    elseif ( ktype == 2 ) then
      dr = sqrt(r2)
      call get_Fab(r, h, fab)
      fab = 0.5*fab

      ! Hes(1,1) = ((dim+2)*r11 - 1)*fab
      ! if ( dim /= 1 ) then
      !   Hes(1,2) = (dim+2)*r12*fab
      !   Hes(2,1) = Hes(1,2)
      !   Hes(2,2) = ((dim+2)*r22 - 1)*fab
      !   if ( dim == 3 ) then
      !     Hes(1,3) = (dim+2)*r13*fab
      !     Hes(3,1) = Hes(1,3)
      !     Hes(2,3) = (dim+2)*r23*fab
      !     Hes(3,2) = Hes(2,3)
      !     Hes(3,3) = ((dim+2)*r33 - 1)*fab
      !   end if
      ! end if
      Hes(1,1) = ((dim+2)*r11 - 1)*fab
      if ( dim /= 1 ) then
        Hes(1,2) = (dim+2)*r12*fab
        Hes(2,1) = Hes(1,2)
        Hes(2,2) = ((dim+2)*r22 - 1)*fab
        if ( dim == 3 ) then
          Hes(1,3) = (dim+2)*r13*fab
          Hes(3,1) = Hes(1,3)
          Hes(2,3) = (dim+2)*r23*fab
          Hes(3,2) = Hes(2,3)
          Hes(3,3) = ((dim+2)*r33 - 1)*fab
        end if
      end if
    elseif ( ktype == 4 ) then
      dr = sqrt(r2)
      call get_FW(r, h, fab)
      fab = 0.5*fab

      Hes(1,1) = ((dim+2)*r11 - 1)*fab
      if ( dim /= 1 ) then
        Hes(1,2) = (dim+2)*r12*fab
        Hes(2,1) = Hes(1,2)
        Hes(2,2) = ((dim+2)*r22 - 1)*fab
        if ( dim == 3 ) then
          Hes(1,3) = (dim+2)*r13*fab
          Hes(3,1) = Hes(1,3)
          Hes(2,3) = (dim+2)*r23*fab
          Hes(3,2) = Hes(2,3)
          Hes(3,3) = ((dim+2)*r33 - 1)*fab
        end if
      end if
    end if
  end subroutine

  pure subroutine get_hessian_rr(r, h, Hes)
  ! subroutine get_hessian_rr(r, h, Hes)
  ! Hessian for artificial viscosity term
  ! H^* = H - F_{ab}/|r_{ab}| I_3       // for momentum methods
  ! H^* = H - f'/ q I_3                 // for direct derivatives
  ! same as original hessian, but without diagonal matrix,
  ! that corresponds to Laplacian term
    real, intent(in)  :: r(3), h
    real, intent(out) :: Hes(3,3)
    real              :: r2, dr, df, ddf, fab, q
    real              :: r11, r12, r13, r22, r23, r33, cstart, dfq
    integer :: ktype

    call getkerntype(ktype)

    r2 = dot_product(r,r)

    r11 = r(1)*r(1)/r2
    r12 = r(1)*r(2)/r2
    r13 = r(1)*r(3)/r2
    r22 = r(2)*r(2)/r2
    r23 = r(2)*r(3)/r2
    r33 = r(3)*r(3)/r2
    cstart = wCv/h**(dim+2)
    Hes(:,:) = 0.

    if (ktype == 1) then
      ! print*, 'Calc with using of dw/dxdy'
      dr = sqrt(r2)
      q = dr / h
      call kddf(q, ddf)
      call kdf(q, df)
      dfq = df/q

      Hes(1,1) = cstart*(ddf - dfq)*r11
      if ( dim /= 1 ) then
        Hes(1,2) = cstart*(ddf - dfq)*r12
        Hes(2,1) = Hes(1,2)
        Hes(2,2) = cstart*(ddf - dfq)*r22
        if ( dim == 3 ) then
          Hes(1,3) = cstart*(ddf - dfq)*r13
          Hes(3,1) = Hes(1,3)
          Hes(2,3) = cstart*(ddf - dfq)*r23
          Hes(3,2) = Hes(2,3)
          Hes(3,3) = cstart*(ddf - dfq)*r33
        end if
      end if
    elseif ( ktype == 2 ) then
      ! print*, 'Calc with using of F_{ab}'
      dr = sqrt(r2)
      call get_Fab(r, h, fab)
      fab = 0.5*fab

      Hes(1,1) = (dim+2)*r11*fab
      if ( dim /= 1 ) then
        Hes(1,2) = (dim+2)*r12*fab
        Hes(2,1) = Hes(1,2)
        Hes(2,2) = (dim+2)*r22*fab
        if ( dim == 3 ) then
          Hes(1,3) = (dim+2)*r13*fab
          Hes(3,1) = Hes(1,3)
          Hes(2,3) = (dim+2)*r23*fab
          Hes(3,2) = Hes(2,3)
          Hes(3,3) = (dim+2)*r33*fab
        end if
      end if
    elseif ( ktype == 4 ) then
      dr = sqrt(r2)
      call get_FW(r, h, fab)
      fab = 0.5*fab

      Hes(1,1) = (dim+2)*r11*fab
      if ( dim /= 1 ) then
        Hes(1,2) = (dim+2)*r12*fab
        Hes(2,1) = Hes(1,2)
        Hes(2,2) = (dim+2)*r22*fab
        if ( dim == 3 ) then
          Hes(1,3) = (dim+2)*r13*fab
          Hes(3,1) = Hes(1,3)
          Hes(2,3) = (dim+2)*r23*fab
          Hes(3,2) = Hes(2,3)
          Hes(3,3) = (dim+2)*r33*fab
        end if
      end if
    end if
  end subroutine
end module
