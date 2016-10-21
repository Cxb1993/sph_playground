module IC
  use kernel
  use BC

  implicit none

  public :: setup

  private
contains

  subroutine setup(t, dim, nx, n, sk, g, pos, vel, acc, mas, den, sln, prs, uie, cf, kcf, dcf)
    character(len=*), intent(in) :: t
    integer, intent(in)  :: nx, dim
    real, intent(in)     :: sk, g
    real, intent(out)    :: pos(nx,3), vel(nx,3), acc(nx,3),&
                            mas(nx), den(nx), sln(nx), prs(nx), uie(nx), cf(nx), kcf(nx), dcf(nx)
    integer, intent(out) :: n
    real                 :: spatVarBrdrs11, spatVarBrdrs12, spatVarBrdrs21, spatVarBrdrs22,&
                            spatVarBrdrs31, spatVarBrdrs32, parSpacing1, parSpacing2,&
                            shockPressure1, shockPressure2, shockDensity1, shockDensity2, &
                            x, y, z, sp
    integer              :: nb, nbnewX1, nbnewY1, nbnewZ1, nbnewX2, nbnewY2, nbnewZ2, &
                            brdarrX1(nx), brdarrY1(nx), brdarrZ1(nx), brdarrX2(nx), &
                            brdarrY2(nx), brdarrZ2(nx)

    call set_dim(dim)
    call set_tasktype(t)

    select case (t)
    case ('hydroshock')
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
    case ('temperhomog01')
      nb = 3
      spatVarBrdrs11 = -1
      spatVarBrdrs12 = 1
      spatVarBrdrs21 = 0.
      spatVarBrdrs22 = 0.
      spatVarBrdrs31 = 0.
      spatVarBrdrs32 = 0.
      if (dim.gt.1) then
        spatVarBrdrs21 = 0
        spatVarBrdrs22 = 0.5
      end if
      if (dim.eq.3) then
        spatVarBrdrs31 = 0.
        spatVarBrdrs32 = 0.05
      end if

      parSpacing1 = 0.025
      parSpacing2 = 0.025

      shockDensity1 = 1000.
      shockDensity2 = 1000.
    end select

    n = 1
    nbnewX1 = 1
    nbnewY1 = 1
    nbnewZ1 = 1
    nbnewX2 = 1
    nbnewY2 = 1
    nbnewZ2 = 1

    x = spatVarBrdrs11
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
          if (x.lt.(spatVarBrdrs11 + nb * sp)) then
            brdarrX1(nbnewX1) = n
            nbnewX1 = nbnewX1 + 1
          else if (x.gt.(spatVarBrdrs12 - nb * sp)) then
            brdarrX2(nbnewX2) = n
            nbnewX2 = nbnewX2 + 1
          end if
          if (dim.gt.1) then
            if (y.lt.(spatVarBrdrs21 + nb * sp)) then
              brdarrY1(nbnewY1) = n
              nbnewY1 = nbnewY1 + 1
            else if (y.gt.(spatVarBrdrs22 - nb * sp)) then
              brdarrY2(nbnewY2) = n
              nbnewY2 = nbnewY2 + 1
            end if
          end if
          select case (t)
          case ('hydroshock')
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
          case ('temperhomog01')
            if (x<0) then
              vel(n,:) = 0.
              acc(n,:) = 0.
              mas(n) = (sp**dim) * shockDensity1
              den(n) = shockDensity1
              sln(n) = sk * sp
              prs(n) = 0
              uie(n) = 0.
              cf(n)  = 0.
              kcf(n) = 1.
            else
              vel(n,:) = 0.
              acc(n,:) = 0.
              mas(n) = (sp**dim) * shockDensity2
              den(n) = shockDensity2
              sln(n) = sk * sp
              prs(n) = 0
              uie(n) = 1.
              cf(n)  = 1.
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
    nbnewX1 = nbnewX2 - 1
    nbnewY1 = nbnewY2 - 1
    nbnewZ1 = nbnewZ2 - 1
    nbnewX2 = nbnewX2 - 1
    nbnewY2 = nbnewY2 - 1
    nbnewZ2 = nbnewZ2 - 1
    n = n - 1

    print *, '#    placed:', n
    print *, '#  border-x:', nbnewX1, nbnewX2
    print *, '#  border-y:', nbnewY1, nbnewY2
    print *, '#  border-z:', nbnewZ1, nbnewZ2

    call set_ns(n)
    call set_border(11, nbnewX1, brdarrX1)
    call set_border(12, nbnewX2, brdarrX2)
    call set_border(21, nbnewY1, brdarrY1)
    call set_border(22, nbnewY2, brdarrY2)
    call set_border(31, nbnewZ1, brdarrZ1)
    call set_border(32, nbnewZ2, brdarrZ2)
  end subroutine setup
end module IC
