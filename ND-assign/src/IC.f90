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
    real, intent(in)     :: pspc1, pspc2, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2
    real, intent(out)    :: pos(3,nx), vel(3,nx), acc(3,nx),&
                            mas(nx), den(nx), sln(nx), prs(nx), uie(nx), cf(nx), kcf(nx), dcf(nx)
    integer, intent(out) :: n
    real                 :: prs1, prs2, rho1, rho2, x, y, z, sp, eps
    integer              :: nb, nbnewX1, nbnewY1, nbnewZ1, nbnewX2, nbnewY2, nbnewZ2, &
                            brdarrX1(nx), brdarrY1(nx), brdarrZ1(nx), brdarrX2(nx), &
                            brdarrY2(nx), brdarrZ2(nx)

    call set_dim(dim)
    call set_tasktype(tt)
    call set_kerntype(kt)
    eps = 1e-8

    nb = 0
    rho1 = 0.
    rho2 = 0.

    select case (tt)
    case ('hydroshock')
      nb = 4
      prs1 = 1.
      prs2 = 0.1
      rho1 = 1.
      rho2 = 0.125
    case ('infslb')
      nb = 1
      rho1 = 1.
      rho2 = 1.
    end select

    n = 1
    nbnewX1 = 1
    nbnewY1 = 1
    nbnewZ1 = 1
    nbnewX2 = 1
    nbnewY2 = 1
    nbnewZ2 = 1

    x = brdx1
    do while ((x >= brdx1).and.(x <= brdx2 + eps))
      if (x.lt.0) then
        sp = pspc1
      else
        sp = pspc2
      end if
      y = brdy1
      do while ((y >= brdy1).and.(y <= brdy2 + eps))
        z = brdz1
        do while ((z >= brdz1).and.(z <= brdz2 + eps))
          pos(1,n) = x
          pos(2,n) = y
          pos(3,n) = z
          if (x.lt.(brdx1 + nb * sp)) then
            brdarrX1(nbnewX1) = n
            nbnewX1 = nbnewX1 + 1
          else if (x.gt.(brdx2 + eps - nb * sp)) then
            brdarrX2(nbnewX2) = n
            nbnewX2 = nbnewX2 + 1
          end if
          if (dim.gt.1) then
            if (y.lt.(brdy1 + nb * sp)) then
              brdarrY1(nbnewY1) = n
              nbnewY1 = nbnewY1 + 1
            else if (y.gt.(brdy2 + eps - nb * sp)) then
              brdarrY2(nbnewY2) = n
              nbnewY2 = nbnewY2 + 1
            end if
            if (dim.eq.3) then
              if (z.lt.(brdz1 + nb * sp)) then
                brdarrZ1(nbnewZ1) = n
                nbnewZ1 = nbnewZ1 + 1
              else if (z.gt.(brdz2 + eps - nb * sp)) then
                brdarrZ2(nbnewZ2) = n
                nbnewZ2 = nbnewZ2 + 1
              end if
            end if
          end if
          !
          ! common values
          !
          vel(:,n) = 0.
          acc(:,n) = 0.
          dcf(n) = 0.
          sln(n) = sk * sp
          prs(n) = 0

          select case (tt)
          case ('hydroshock')
            if (x<0) then
              mas(n) = (sp**dim) * rho1
              den(n) = rho1
              prs(n) = prs1
              uie(n) = prs1 / (g - 1) / rho1
            else
              mas(n) = (sp**dim) * rho2
              den(n) = rho2
              prs(n) = prs2
              uie(n) = prs2 / (g - 1) / rho2
            end if
          case ('infslb')
            if (x<0) then
              mas(n) = (sp**dim) * rho1
              den(n) = rho1
              ! cf(n)  = sin(pi * (x - brdx1) / abs(brdx2-brdx1))*sin(pi * (y - brdy1) / abs(brdy2-brdy1))
              ! cf(n)  = sin(2 * pi * (x + 1.) / abs(brdx2-brdx1))
              cf(n)  = 0.
              kcf(n) = 1.
            else if (x > 0+eps) then
              mas(n) = (sp**dim) * rho2
              den(n) = rho2
              prs(n) = 0
              ! cf(n)  = sin(pi * (x - brdx1) / abs(brdx2-brdx1))*sin(pi * (y - brdy1) / abs(brdy2-brdy1))
              ! cf(n)  = sin(2 * pi * (x + 1.) / abs(brdx2-brdx1))
              cf(n)  = 1.
              kcf(n) = 1.
            else
              mas(n) = (sp**dim) * rho2
              den(n) = rho2
              ! cf(n)  = sin(pi * (x - brdx1) / abs(brdx2-brdx1))*sin(pi * (y - brdy1) / abs(brdy2-brdy1))
              ! cf(n)  = sin(pi * (x + 1.) / abs(brdx2-brdx1))
              cf(n)  = .5
              kcf(n) = 1.
            end if
            uie(n) = cf(n) / cv
          end select
          z = z + sp
          n = n + 1
        end do
        y = y + sp
      end do
      x = x + sp
    end do
    nbnewX1 = nbnewX1 - 1
    nbnewY1 = nbnewY1 - 1
    nbnewZ1 = nbnewZ1 - 1
    nbnewX2 = nbnewX2 - 1
    nbnewY2 = nbnewY2 - 1
    nbnewZ2 = nbnewZ2 - 1
    n = n - 1

    print *, '#    placed:', n
    print *, '#  border-x:', nbnewX1, nbnewX2
    print *, '#  border-y:', nbnewY1, nbnewY2
    print *, '#  border-z:', nbnewZ1, nbnewZ2

    call set_particles_numbers(n, nb)
    call set_border(11, nbnewX1, brdarrX1)
    call set_border(12, nbnewX2, brdarrX2)
    call set_border(21, nbnewY1, brdarrY1)
    call set_border(22, nbnewY2, brdarrY2)
    call set_border(31, nbnewZ1, brdarrZ1)
    call set_border(32, nbnewZ2, brdarrZ2)
  end subroutine setup
end module IC
