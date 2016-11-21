module IC
  use kernel
  use BC

  implicit none

  public :: setup

  private
  real, parameter     :: pi = 4.*atan(1.)

contains

  subroutine setup(tt, kt, dim, nx, n, sk, g, cv, &
    pspc1, pspc2, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, &
    pos, vel, acc, mas, den, sln, prs, uie, cf, kcf, dcf)
    character(len=*), intent(in) :: tt, kt
    integer, intent(in)  :: nx, dim
    real, intent(in)     :: sk, g, cv
    real, intent(in)     :: brdx1, brdx2, brdy1, brdy2, brdz1, brdz2
    real, intent(out)    :: pos(3,nx), vel(3,nx), acc(3,nx),&
                            mas(nx), den(nx), sln(nx), prs(nx), uie(nx), cf(nx), kcf(nx), dcf(nx)
    real, intent(inout)  :: pspc1, pspc2
    integer, intent(out) :: n
    real                 :: prs1, prs2, rho1, rho2, k1, k2, x, y, z, sp, spx, spy, spz
    integer              :: i, j, k, nb, ix, iy, iz, bdx, bdy, bdz, ibx, iby, ibz, &
                            nbnewX1, nbnewY1, nbnewZ1, nbnewX2,&
                            nbnewY2, nbnewZ2, brdarrX1(nx), brdarrY1(nx), brdarrZ1(nx), brdarrX2(nx), &
                            brdarrY2(nx), brdarrZ2(nx)

    call set_dim(dim)
    call set_tasktype(tt)
    call set_kerntype(kt)

    nb = 0
    rho1 = 1.
    rho2 = 1.
    prs1 = 1.
    prs2 = 1.
    k1 = 1.
    k2 = 1.

    select case (tt)
    case ('hydroshock')
      nb = 4
      prs1 = 1.
      prs2 = 0.1
      rho1 = 1.
      rho2 = 0.125
    case ('infslb')
      nb = 1
    case ('hc-sinx')
      nb = 3
    end select

    n = 1
    nbnewX1 = 1
    nbnewY1 = 1
    nbnewZ1 = 1
    nbnewX2 = 1
    nbnewY2 = 1
    nbnewZ2 = 1

    if (pspc1 /= pspc2) then
      ! if (dim > 1) then
      !   nby = nb
      !   if (dim == 3) then
      !     nbz = nb
      !   end if
      ! end if
      x = brdx1! - pspc1 * nbx
      ! do while ((x >= brdx1 - pspc1 * nbx).and.(x <= brdx2 + pspc2 * nbx))
      do while ((x >= brdx1).and.(x <= brdx2))
        if (x < 0) then
          sp = pspc1
        else
          sp = pspc2
        end if
        y = brdy1! - pspc1 * nby
        ! do while ((y >= brdy1 - pspc1 * nby).and.(y <= brdy2 + pspc2 * nby))
        do while ((y >= brdy1).and.(y <= brdy2))
          z = brdz1! - pspc1 * nbz
          ! do while ((z >= brdz1 - pspc1 * nbz).and.(z <= brdz2 + pspc2 * nbz))
          do while ((z >= brdz1).and.(z <= brdz2))
            pos(1,n) = x
            pos(2,n) = y
            pos(3,n) = z
            if (x < brdx1 + nb * sp) then
              brdarrX1(nbnewX1) = n
              nbnewX1 = nbnewX1 + 1
            else if (x > brdx2 - nb * sp) then
              brdarrX2(nbnewX2) = n
              nbnewX2 = nbnewX2 + 1
            end if
            if (dim > 1) then
              ! if (y < brdy1) then
              if (y < brdy1 + nb * sp) then
                brdarrY1(nbnewY1) = n
                nbnewY1 = nbnewY1 + 1
              ! else if (y > brdy2) then
              else if (y > brdy2 - nb * sp) then
                brdarrY2(nbnewY2) = n
                nbnewY2 = nbnewY2 + 1
              end if
              if (dim == 3) then
                if (z < brdz1 + nb * sp) then
                  brdarrZ1(nbnewZ1) = n
                  nbnewZ1 = nbnewZ1 + 1
                else if (z > brdz2 - nb * sp) then
                  brdarrZ2(nbnewZ2) = n
                  nbnewZ2 = nbnewZ2 + 1
                end if
              end if
            end if
            z = z + sp
            n = n + 1
          end do
          y = y + sp
        end do
        x = x + sp
      end do
    else
      ix = int((brdx2-brdx1)/pspc1)
      spx = merge(0.,(brdx2-brdx1)/ix, ix == 0)
      pspc1 = spx
      pspc2 = spx
      iy = int((brdy2-brdy1)/pspc1)
      spy = merge(0.,(brdy2-brdy1)/iy, iy == 0)
      iz = int((brdz2-brdz1)/pspc1)
      spz = merge(0.,(brdz2-brdz1)/iz, iz == 0)
      if (nb > 0) then
        bdx = nb
        bdy = merge(nb, 0, dim > 1)
        bdz = merge(nb, 0, dim == 3)
        ibx = 0
        iby = 0
        ibz = 0
      else
        bdx = 0
        bdy = 0
        bdz = 0
        ibx = abs(nb)
        iby = abs(nb)
        ibz = abs(nb)
      end if
      do i = (0-bdx),(ix+bdx)
        x = brdx1 + i * spx
        do j = (0-bdy),(iy+bdy)
          y = brdy1 + j * spy
          do k = (0-bdz),(iz+bdz)
            z = brdz1 + k * spz
            if (i < 0 + ibx) then
              brdarrX1(nbnewX1) = n
              nbnewX1 = nbnewX1 + 1
            else if (i > ix - ibx) then
              brdarrX2(nbnewX2) = n
              nbnewX2 = nbnewX2 + 1
            end if
            if (dim > 1) then
              if (j < 0 + iby) then
                brdarrY1(nbnewY1) = n
                nbnewY1 = nbnewY1 + 1
              else if (j > iy - iby) then
                brdarrY2(nbnewY2) = n
                nbnewY2 = nbnewY2 + 1
              end if
              if (dim == 3) then
                if (k < 0 + ibz) then
                  brdarrZ1(nbnewZ1) = n
                  nbnewZ1 = nbnewZ1 + 1
                else if (k > iz - ibz) then
                  brdarrZ2(nbnewZ2) = n
                  nbnewZ2 = nbnewZ2 + 1
                end if
              end if
            end if
            pos(1,n) = x
            pos(2,n) = y
            pos(3,n) = z
            n = n + 1
          end do
        end do
      end do
    end if

    nbnewX1 = nbnewX1 - 1
    nbnewY1 = nbnewY1 - 1
    nbnewZ1 = nbnewZ1 - 1
    nbnewX2 = nbnewX2 - 1
    nbnewY2 = nbnewY2 - 1
    nbnewZ2 = nbnewZ2 - 1
    n = n - 1

    write(*, "(A, F7.5, A, F7.5)") " # actual dx:   x1=", pspc1, "   x2=", pspc2
    print *, '#    placed:', n
    print *, '#      real:', n - nbnewX1 - nbnewX2 - nbnewY1 - nbnewY2 - nbnewZ1 - nbnewZ2
    print *, '#  border-x:', nbnewX1, nbnewX2
    print *, '#  border-y:', nbnewY1, nbnewY2
    print *, '#  border-z:', nbnewZ1, nbnewZ2

    call set_particles_numbers(n, abs(nb))
    call set_border(11, nbnewX1, brdarrX1)
    call set_border(12, nbnewX2, brdarrX2)
    call set_border(21, nbnewY1, brdarrY1)
    call set_border(22, nbnewY2, brdarrY2)
    call set_border(31, nbnewZ1, brdarrZ1)
    call set_border(32, nbnewZ2, brdarrZ2)
    !
    ! common values
    !
    do i=1,n
      sp = merge(pspc1, pspc2, pos(1,i) < 0)

      vel(:,i) = 0.
      acc(:,i) = 0.
      dcf(i) = 0.
      sln(i) = sk * sp
      prs(i) = 0

      if (pos(1,i) < 0) then
        mas(i) = (sp**dim) * rho1
        den(i) = rho1
        prs(i) = prs1
        kcf(i) = k1
      else
        mas(i) = (sp**dim) * rho2
        den(i) = rho2
        prs(i) = prs2
        kcf(i) = k2
      end if

      select case (tt)
      case ('hydroshock')
        uie(i) = merge(prs1/(g-1)/rho1, prs2/(g-1)/rho2, pos(1,i) < 0)
      case ('infslb')
        cf(i) = merge(0., 1., pos(1,i) < 0)
        uie(i) = cf(i) / cv
      case ('hc-sinx')
        cf(i)  = sin(2 * pi * (pos(1,i) + 1.) / abs(brdx2-brdx1))
        uie(i) = cf(i) / cv
      end select
    end do
  end subroutine setup
end module IC
