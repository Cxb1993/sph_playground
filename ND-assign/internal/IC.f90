module IC
  use kernel
  use BC

  implicit none

  public :: setup

  private
  real, parameter     :: pi = 4.*atan(1.)

contains

  subroutine setup(t, dim, nx, n, sk, g, cv, pos, vel, acc, mas, den, sln, prs, uie, cf, kcf, dcf)
    character(len=*), intent(in) :: t
    integer, intent(in)  :: nx, dim
    real, intent(in)     :: sk, g, cv
    real, intent(out)    :: pos(nx,3), vel(nx,3), acc(nx,3),&
                            mas(nx), den(nx), sln(nx), prs(nx), uie(nx), cf(nx), kcf(nx), dcf(nx)
    integer, intent(out) :: n
    real                 :: brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, parSpacing1, parSpacing2,&
                            shockPressure1, shockPressure2, shockDensity1, shockDensity2, &
                            x, y, z, sp, eps
    integer              :: nb, nbnewX1, nbnewY1, nbnewZ1, nbnewX2, nbnewY2, nbnewZ2, &
                            brdarrX1(nx), brdarrY1(nx), brdarrZ1(nx), brdarrX2(nx), &
                            brdarrY2(nx), brdarrZ2(nx)

    call set_dim(dim)
    call set_tasktype(t)
    eps = 0.00001

    nb = 0
    brdx1 = 0.
    brdx2 = 0.
    brdy1 = 0.
    brdy2 = 0.
    brdz1 = 0.
    brdz2 = 0.
    shockDensity1 = 0.
    shockDensity2 = 0.
    parSpacing1 = 0.
    parSpacing2 = 0.

    select case (t)
    case ('hydroshock')
      nb = 4
      brdx1 = -0.5
      brdx2 = 0.5
      brdy1 = 0.
      brdy2 = 0.
      brdz1 = 0.
      brdz2 = 0.
      if (dim.gt.1) then
        brdy1 = 0.
        brdy2 = 0.5
        if (dim.eq.3) then
          brdz1 = 0.
          brdz2 = 0.05
        end if
      end if
      parSpacing1 = 0.001
      parSpacing2 = 0.008
      shockPressure1 = 1.
      shockPressure2 = 0.1
      shockDensity1 = 1.
      shockDensity2 = 0.125
    case ('heatslab')
      nb = 1
      brdx1 = -1.
      brdx2 = 1.
      brdy1 = 0.
      brdy2 = 0.
      brdz1 = 0.
      brdz2 = 0.
      if (dim.gt.1) then
        brdy1 = -1.
        brdy2 = 1.
        if (dim.eq.3) then
          brdz1 = 0.
          brdz2 = 0.05
        end if
      end if
      parSpacing1 = 0.014
      parSpacing2 = 0.014
      shockDensity1 = 1.
      shockDensity2 = 1.
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
        sp = parSpacing1
      else
        sp = parSpacing2
      end if
      y = brdy1
      do while ((y >= brdy1).and.(y <= brdy2 + eps))
        z = brdz1
        do while ((z >= brdz1).and.(z <= brdz2 + eps))
          pos(n,1) = x
          pos(n,2) = y
          pos(n,3) = z
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
          end if
          select case (t)
          case ('hydroshock')
            if (x<0) then
              vel(n,:) = 0.
              acc(n,:) = 0.
              dcf(n) = 0.
              mas(n) = (sp**dim) * shockDensity1
              den(n) = shockDensity1
              sln(n) = sk * sp
              prs(n) = shockPressure1
              uie(n) = shockPressure1 / (g - 1) / shockDensity1
            else
              vel(n,:) = 0.
              acc(n,:) = 0.
              dcf(n) = 0.
              mas(n) = (sp**dim) * shockDensity2
              den(n) = shockDensity2
              sln(n) = sk * sp
              prs(n) = shockPressure2
              uie(n) = shockPressure2 / (g - 1) / shockDensity2
            end if
          case ('heatslab')
            if (x<0) then
              vel(n,:) = 0.
              acc(n,:) = 0.
              dcf(n) = 0.
              mas(n) = (sp**dim) * shockDensity1
              den(n) = shockDensity1
              sln(n) = sk * sp
              prs(n) = 0
              ! cf(n)  = sin(pi * (x - brdx1) / abs(brdx2-brdx1))*sin(pi * (y - brdy1) / abs(brdy2-brdy1))
              ! cf(n)  = sin(pi * (x + 1.) / abs(brdx2-brdx1))
              cf(n)  = 1.
              uie(n) = cf(n) / cv
              kcf(n) = 1000.
            else
              vel(n,:) = 0.
              acc(n,:) = 0.
              dcf(n) = 0.
              mas(n) = (sp**dim) * shockDensity2
              den(n) = shockDensity2
              sln(n) = sk * sp
              prs(n) = 0
              ! cf(n)  = sin(pi * (x - brdx1) / abs(brdx2-brdx1))*sin(pi * (y - brdy1) / abs(brdy2-brdy1))
              ! cf(n)  = sin(pi * (x + 1.) / abs(brdx2-brdx1))
              cf(n)  = 2.
              uie(n) = cf(n) / cv
              kcf(n) = 1.
            end if
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
