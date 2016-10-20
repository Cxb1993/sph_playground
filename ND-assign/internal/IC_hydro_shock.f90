module IC_hydro_shock
  use kernel
  use BC

  implicit none

  public :: setup_hydro_shock

  private
contains

  subroutine setup_hydro_shock(dim, nx, n, sk, g, pos, vel, acc, mas, den, sln, prs, uie)
    integer, intent(in)  :: nx, dim
    real, intent(in)     :: sk, g
    real, intent(out)    :: pos(nx,3), vel(nx,3), acc(nx,3), mas(nx), den(nx), sln(nx), prs(nx), uie(nx)
    integer, intent(out) :: n
    real                 :: spatVarBrdrs11, spatVarBrdrs12, spatVarBrdrs21, spatVarBrdrs22, spatVarBrdrs31, spatVarBrdrs32
    real                 :: parSpacing1, parSpacing2, shockPressure1, shockPressure2, shockDensity1, shockDensity2
    real                 :: x, y, z, sp
    integer              :: nb, nbnewX, nbnewY, nbnewZ, brdarrX(nx), brdarrY(nx), brdarrZ(nx)

    call set_dim(dim)

    nb = 3
    spatVarBrdrs11 = -0.5
    spatVarBrdrs12 = 0.5
    spatVarBrdrs21 = 0.
    spatVarBrdrs22 = 0.
    spatVarBrdrs31 = 0.
    spatVarBrdrs32 = 0.
    if (dim.gt.1) then
      spatVarBrdrs21 = 0.
      spatVarBrdrs22 = 0.05
    end if
    if (dim.eq.3) then
      spatVarBrdrs31 = 0.
      spatVarBrdrs32 = 0.05
    end if

    parSpacing1 = 0.001
    parSpacing2 = 0.008

    shockPressure1 = 1.
    shockPressure2 = 0.1

    shockDensity1 = 1.
    shockDensity2 = 0.125

    x = spatVarBrdrs11
    n = 1
    nbnewX = 1
    nbnewY = 1
    nbnewZ = 1
    do while ((x >= spatVarBrdrs11).and.(x <= spatVarBrdrs12))
      if (x.lt.0) then
        sp = parSpacing1
      else
        sp = parSpacing2
      end if
      y = spatVarBrdrs21
      do while ((y >= spatVarBrdrs21).and.(y <= spatVarBrdrs22))
        z = spatVarBrdrs31
        do while ((z >= spatVarBrdrs31).and.(z <= spatVarBrdrs32))
          pos(n,1) = x
          pos(n,2) = y
          pos(n,3) = z
          if ((x.lt.(spatVarBrdrs11 + nb * sp)).or.(x.gt.(spatVarBrdrs12 - nb * sp))) then
            brdarrX(nbnewX) = n
            nbnewX = nbnewX + 1
          end if
          if (dim.gt.1) then
            if ((y.lt.(spatVarBrdrs21 + nb * sp)).or.(y.gt.(spatVarBrdrs22 - nb * sp))) then
              brdarrY(nbnewY) = n
              nbnewY = nbnewY + 1
            end if
          end if

          if (x<0) then
            vel(n,:) = 0.
            acc(n,:) = 0.
            mas(n) = (sp**dim) * shockDensity1
            den(n) = shockDensity1
            sln(n) = sk * sp
            prs(n) = shockPressure1
            uie(n) = shockPressure1 / (g - 1) / shockDensity1
          else
            vel(n,:) = 0.
            acc(n,:) = 0.
            mas(n) = (sp**dim) * shockDensity2
            den(n) = shockDensity2
            sln(n) = sk * sp
            prs(n) = shockPressure2
            uie(n) = shockPressure2 / (g - 1) / shockDensity2
          end if
          z = z + sp
          n = n + 1
        end do
        y = y + sp
      end do
      x = x + sp
    end do
    nbnewX = nbnewX - 1
    nbnewY = nbnewY - 1
    nbnewZ = nbnewZ - 1
    n = n - 1

    print *, '#    placed:', n
    print *, '#  border-x:', nbnewX
    print *, '#  border-y:', nbnewY
    print *, '#  border-z:', nbnewZ

    call set_ns(n)
    call set_borders(nbnewX, nbnewY, nbnewZ,&
                     brdarrX(1:nbnewX), brdarrY(1:nbnewY), brdarrZ(1:nbnewZ))
  end subroutine setup_hydro_shock
end module IC_hydro_shock
