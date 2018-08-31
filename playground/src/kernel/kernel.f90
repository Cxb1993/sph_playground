module kernel
  use errprinter, only: error
  use const
  use kernel_base
  use state,      only: gcoordsys, &
                        getddwtype, &
                        setdim

  implicit none

  public :: nw, get_dw_dh, get_w, initkernel, &
            n2w, get_krad, hessian,&
            hessian_rr, getneibnumber, precalcKernel
  save

  integer :: dim

  procedure (ainw),     pointer :: nw         => null()
  procedure (ain2w),    pointer :: n2w        => null()
  procedure (aihesw),   pointer :: hessian    => null()
  procedure (aihesrrw), pointer :: hessian_rr => null()

  abstract interface
    subroutine ainw(rab, ra, rb, dr, h, onw)
      real, intent(out) :: onw(3)
      real, intent (in) :: rab(3), ra(3), rb(3), dr, h
    end subroutine ainw
  end interface

  abstract interface
    subroutine aihesw(rab, ra, rb, h, hesw)
      real, intent(out) :: hesw(3,3)
      real, intent (in) :: rab(3), ra(3), rb(3), h
    end subroutine aihesw
  end interface

  abstract interface
    pure subroutine aihesrrw(rab, h, hesrrw)
      real, intent(out) :: hesrrw(3,3)
      real, intent (in) :: rab(3), h
    end subroutine aihesrrw
  end interface

  abstract interface
    pure subroutine ain2w(rab, h, on2w)
      real, intent(out) :: on2w
      real, intent (in) :: rab(3), h
    end subroutine ain2w
  end interface

  private
  integer, parameter :: kernelres=1001
  real :: sf(kernelres), sdf(kernelres), sddf(kernelres), dk

contains
  pure subroutine get_krad(kr)
    real, intent(out) :: kr
    kr = krad
  end subroutine get_krad

  subroutine getneibnumber(nn)
    integer, intent(out) :: nn

    nn = returnneibnum
  end subroutine

  subroutine initkernel(indim)
    integer, intent(in) :: indim
    integer :: cs, kt

    dim = indim
    call setdimbase(dim)
    call setdim(dim)

    call gcoordsys(cs)
    call getddwtype(kt)
    ! call precalcKernel()

    if (cs == 1) then
      nw => nw_cart
      hessian_rr => hessian_rr_fab_cart
      ! nw => nw_cart_precalc
      if (kt == esd_n2w) then
        hessian => hessian_ddw_cart
        n2w     => n2w_cart
      else if (kt == esd_fab) then
        hessian => hessian_fab_cart
        n2w     => Fab_cart
      else if ((kt == esd_2nw_ds).or.(kt == esd_2nw_sd)) then
        hessian => null()
        n2w     => null()
      else if (kt == esd_fw) then
        hessian => hessian_fw_cart
        n2w     => FW_cart
      else
        error stop "Wrong kernel type in kernel init."
      end if
    else if (cs == 2) then
      nw => nw_cyl

      if (kt == esd_n2w) then
        hessian => hessian_ddw_cyl
      else if (kt == esd_fab) then
        ! hessian => hessian_fab_cart
      else if ((kt == esd_2nw_ds).or.(kt == esd_2nw_sd)) then
        ! no more actions
      else if (kt == esd_fw) then
        ! hessian => hessian_fw_cart
      else
        call error('Wrong kernel type', '', __FILE__, __LINE__)
      end if
    else
      call error('Wrong coordinate system', '', __FILE__, __LINE__)
    end if
  end subroutine

 subroutine precalcKernel()
   real    :: q, dq, f
   integer :: i

   dq = krad/(kernelres-1)
   dk = dq
   do i = 1,kernelres
     q = dq*(i-1)
     call kf(q, f)
     sf(i) = wCv * f
     call kdf(q, f)
     sdf(i) = wCv * f
     call kddf(q, f)
     sddf(i) = wCv * f
     ! print*, i, nint(q/dk)+1, dq*(i-1), sf(i), sdf(i), sddf(i)
   end do
   ! read*
 end subroutine precalcKernel

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

  pure subroutine nw_cart(rab, ra, rb, dr, h, nw)
    real, intent(in)  :: rab(3), ra(3), rb(3), dr, h
    real, intent(out) :: nw(3)
    real              :: df, q

    q = dr / h
    call kdf(q, df)

    nw(:) = wCv * df * rab(:) / h**(dim+2) / q
  end subroutine nw_cart

  pure subroutine nw_cart_precalc(rab, ra, rb, dr, h, nw)
    real, intent(in)  :: rab(3), ra(3), rb(3), dr, h
    real, intent(out) :: nw(3)
    real              :: q

    q = dr / h
    nw(:) = sdf(nint(q/dk)+1) * rab(:) / h**(dim+2) / q
  end subroutine nw_cart_precalc

  pure subroutine nw_cyl(rab, ra, rb, dr, h, nw)
    ! get cylindrical nabla kernel for cartesian input
    real, intent(in)  :: rab(3), ra(3), rb(3), dr, h
    real, intent(out) :: nw(3)
    real              :: df, q, cca(3), ccb(3)

    cca(1) = sqrt(ra(1)*ra(1) + ra(2)*ra(2))
    cca(2) = atan(ra(2),ra(1))
    cca(3) = ra(3)

    ccb(1) = sqrt(rb(1)*rb(1) + rb(2)*rb(2))
    ccb(2) = atan(rb(2),rb(1))
    ccb(3) = rb(3)

    q = dr / h
    call kdf(q, df)

    nw(1) = wCv / h**(dim+2) / q * df * (cca(1)-ccb(1)+ccb(1)*(1-cos(cca(2)-ccb(2))))
    if (cca(1) > eps0) then
      nw(2) = wCv / h**(dim+2) / q * df * (cca(1)*ccb(1)*sin(cca(2)-ccb(2))) / cca(1)
    else
      nw(2) = 0.
    end if
    nw(3) = wCv / h**(dim+2) / q * df * (cca(3)-ccb(3))
  end subroutine nw_cyl

  pure subroutine get_dw_dh(r, h, dwdh)
    real, intent(in)  :: r, h
    real, intent(out) :: dwdh
    real              :: f, df, q
    q = r / h

    call kf(q, f)
    call kdf(q, df)

    dwdh = - wCv / h**(dim + 1) * (dim * f + q * df)
  end subroutine get_dw_dh

  pure subroutine FW_cart(r, h, fw)
    real, intent(in)  :: r(3), h
    real, intent(out) :: fw
    real              :: f, dr, q

    dr = sqrt(dot_product(r,r))
    q = dr / h
    call kf(q, f)
    fw = fwc * wCv * f / h ** (dim + 2)
  end subroutine FW_cart

  pure subroutine Fab_cart(r, h, Fab)
    real, intent(in)  :: r(3), h
    real, intent(out) :: Fab
    real              :: q, df

    q = sqrt(dot_product(r,r)) / h
    call kdf(q, df)
    Fab = -2. * wCv * df / h**(dim+2) / q
  end subroutine Fab_cart

  pure subroutine n2w_cart(r, h, n2w)
    real, intent(in)  :: r(3), h
    real, intent(out) :: n2w
    real              :: df, ddf, q

    q = sqrt(dot_product(r,r)) / h

    call kddf(q, ddf)
    call kdf(q, df)
    n2w = wCv*(ddf + (dim - 1) * df / q)/h**(dim + 2)
  end subroutine n2w_cart

  pure subroutine hessian_rr_fab_cart(r, h, Hes)
    ! subroutine get_hessian_rr(r, h, Hes)
    ! Hessian for artificial viscosity term
    ! H^* = H - F_{ab}/|r_{ab}| I_3       // for momentum methods
    ! H^* = H - f'/ q I_3                 // for direct derivatives
    ! same as original hessian, but without diagonal matrix,
    ! that corresponds to Laplacian term
    real, intent(in)  :: r(3), h
    real, intent(out) :: Hes(3,3)
    real              :: r2, dr, fab
    real              :: r11, r12, r13, r22, r23, r33, cstart

    r2 = dot_product(r,r)

    r11 = r(1)*r(1)/r2
    r12 = r(1)*r(2)/r2
    r13 = r(1)*r(3)/r2
    r22 = r(2)*r(2)/r2
    r23 = r(2)*r(3)/r2
    r33 = r(3)*r(3)/r2
    cstart = wCv/h**(dim+2)
    Hes(:,:) = 0.

    dr = sqrt(r2)
    call Fab_cart(r, h, fab)
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
  end subroutine hessian_rr_fab_cart

  pure subroutine hessian_rr_n2w_cart(r, h, Hes)
    ! subroutine get_hessian_rr(r, h, Hes)
    ! Hessian for artificial viscosity term
    ! H^* = H - F_{ab}/|r_{ab}| I_3       // for brookshaw methods
    ! H^* = H - f'/ q I_3                 // for direct derivatives
    ! same as original hessian, but without diagonal matrix,
    ! that corresponds to Laplacian term
    real, intent(out) :: Hes(3,3)
    real, intent(in)  :: r(3), h
    real              :: r2, dr, df, ddf, q
    real              :: r11, r12, r13, r22, r23, r33, cstart, dfq

    r2 = dot_product(r,r)

    r11 = r(1)*r(1)/r2
    r12 = r(1)*r(2)/r2
    r13 = r(1)*r(3)/r2
    r22 = r(2)*r(2)/r2
    r23 = r(2)*r(3)/r2
    r33 = r(3)*r(3)/r2
    cstart = wCv/h**(dim+2)
    Hes(:,:) = 0.

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
  end subroutine hessian_rr_n2w_cart

  subroutine hessian_ddw_cart(rab, ra, rb, h, Hes)
  ! subroutine hessian_ddw_cart(rab, ra, rb, h, Hes)
    real, intent(in)  :: rab(3), ra(3), rb(3), h
    real, intent(out) :: Hes(3,3)
    real              :: r2, dr, df, ddf, q
    real              :: r11, r12, r13, r22, r23, r33, cstart, dfq

    r2 = dot_product(rab,ra(:)-rb(:))

    r11 = rab(1)*rab(1)/r2
    r12 = rab(1)*rab(2)/r2
    r13 = rab(1)*rab(3)/r2
    r22 = rab(2)*rab(2)/r2
    r23 = rab(2)*rab(3)/r2
    r33 = rab(3)*rab(3)/r2
    cstart = wCv/h**(dim+2)
    Hes(:,:) = 0.

    dr = sqrt(r2)

    q = dr / h

    call kddf(q, ddf)
    call kdf(q, df)

    dfq = df/q

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
  end subroutine

  subroutine hessian_w_cart(rab, ra, rb, h, Hes)
    real, intent(in)  :: rab(3), ra(3), rb(3), h
    real, intent(out) :: Hes(3,3)
    real              :: r2, dr, df, ddf, q
    real              :: r11, r12, r13, r22, r23, r33, cstart, dfq

    r2 = dot_product(rab,ra(:)-rb(:))

    r11 = rab(1)*rab(1)/r2
    r12 = rab(1)*rab(2)/r2
    r13 = rab(1)*rab(3)/r2
    r22 = rab(2)*rab(2)/r2
    r23 = rab(2)*rab(3)/r2
    r33 = rab(3)*rab(3)/r2
    cstart = wCv/h**(dim+2)
    Hes(:,:) = 0.

    dr = sqrt(r2)

    q = dr / h

    call kddf(q, ddf)
    call kdf(q, df)

    dfq = df/q

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
  end subroutine

  pure subroutine hessian_fab_cart(rab, ra, rb, h, Hes)
    real, intent(in)  :: rab(3), ra(3), rb(3), h
    real, intent(out) :: Hes(3,3)
    real              :: r2, dr, fab
    real              :: r11, r12, r13, r22, r23, r33, cstart!, a

    r2 = dot_product(rab,ra(:)-rb(:))

    r11 = rab(1)*rab(1)/r2
    r12 = rab(1)*rab(2)/r2
    r13 = rab(1)*rab(3)/r2
    r22 = rab(2)*rab(2)/r2
    r23 = rab(2)*rab(3)/r2
    r33 = rab(3)*rab(3)/r2
    cstart = wCv/h**(dim+2)
    Hes(:,:) = 0.

    dr = sqrt(r2)
    call Fab_cart(rab, h, fab)
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
    ! a = 2./5.
    ! Hes(:,:) = a*Hes(:,:)
    ! Hes(1,1) = Hes(1,1) + 1./3.*(1 - a)
    ! if (dim > 1) then
    !   Hes(2,2) = Hes(2,2) + 1./3.*(1 - a)
    !   if (dim == 3) then
    !     Hes(3,3) = Hes(3,3) + 1./3.*(1 - a)
    !   end if
    ! end if
  end subroutine

  pure subroutine hessian_fw_cart(rab, ra, rb, h, Hes)
    real, intent(in)  :: rab(3), ra(3), rb(3), h
    real, intent(out) :: Hes(3,3)
    real              :: r2, dr, fab
    real              :: r11, r12, r13, r22, r23, r33, cstart

    r2 = dot_product(rab,ra(:) - rb(:))

    r11 = rab(1)*rab(1)/r2
    r12 = rab(1)*rab(2)/r2
    r13 = rab(1)*rab(3)/r2
    r22 = rab(2)*rab(2)/r2
    r23 = rab(2)*rab(3)/r2
    r33 = rab(3)*rab(3)/r2
    cstart = wCv/h**(dim+2)
    Hes(:,:) = 0.

    dr = sqrt(r2)
    call FW_cart(rab, h, fab)
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
  end subroutine

  pure subroutine hessian_ddw_cyl(rab, ra, rb, h, Hes)
    real, intent(in)  :: rab(3), ra(3), rb(3), h
    real, intent(out) :: Hes(3,3)
    real              :: r2, ir2, dr, df, ddf, q, nw(3)
    real              :: r11, r12, r13, r22, r23, r33, cstart, dfq
    real              :: d11, d12, d13, d22, d23, d33
    real              :: cosdphi, sindphi, term1
    real              :: cca(3), ccb(3)

    cca(1) = sqrt(ra(1)*ra(1) + ra(2)*ra(2))
    cca(2) = atan(ra(2),ra(1))
    cca(3) = ra(3)

    ccb(1) = sqrt(rb(1)*rb(1) + rb(2)*rb(2))
    ccb(2) = atan(rb(2),rb(1))
    ccb(3) = rb(3)

    r2 = dot_product(rab,rab)
    cstart = wCv/h**(dim+2)
    dr = sqrt(r2)
    ir2 = 1./r2
    q = dr / h

    cosdphi = cos(cca(2)-ccb(2))
    sindphi = sin(cca(2)-ccb(2))
    term1 = ((cca(1) - ccb(1)) + ccb(1)*(1 - cosdphi))

    r11 = ir2*term1*term1
    d11 = 1.
    r12 = ir2*sindphi*cca(1)*ccb(1)*term1
    d12 = sindphi
    r13 = rab(3)*ir2*term1
    d13 = 0.
    r22 = ir2*cca(1)*cca(1)*ccb(1)*ccb(1)*sindphi*sindphi
    d22 = cca(1)*ccb(1)*cosdphi/q
    r23 = cca(1)*ccb(1)*sindphi*rab(3)*ir2
    d23 = 0.
    r33 = rab(3)*rab(3)*ir2
    d33 = 1.

    Hes(:,:) = 0.

    call kddf(q, ddf)
    call kdf(q, df)
    call nw_cyl(rab, ra, rb, dr, h, nw)

    dfq = df/q
    if (cca(1) > eps0) then
      Hes(1,1) = cstart*((ddf - dfq)*r11 + dfq*d11) + nw(1) / cca(1)
      if (dim /= 1) then
        Hes(1,2) = cstart*((ddf - dfq)*r12 + dfq*d12) / cca(1)
        Hes(2,1) = Hes(1,2)
        Hes(2,2) = cstart*((ddf - dfq)*r22 + dfq*d22) / cca(1)
        if ( dim == 3 ) then
          Hes(1,3) = cstart*((ddf - dfq)*r13 + dfq*d13)
          Hes(3,1) = Hes(1,3)
          Hes(2,3) = cstart*((ddf - dfq)*r23 + dfq*d23) / cca(1)
          Hes(3,2) = Hes(2,3)
          Hes(3,3) = cstart*((ddf - dfq)*r33 + dfq*d33)
        end if
      end if
    end if
  end subroutine hessian_ddw_cyl
end module
