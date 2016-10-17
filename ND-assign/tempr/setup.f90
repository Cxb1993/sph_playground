module tempr_setup
  use kernel

  implicit none

  public :: tempr_homog01, tempr_set_fixed1, tempr_set_fixed3

  private
    integer, save        :: ns
    integer, allocatable :: borderX(:), borderY(:), borderZ(:)

contains

  subroutine tempr_homog01(dim, nx, n, sk, g, pos, vel, acc, mas, den, sln, prs, uie, cf, kcf)
    integer, intent(in)  :: nx, dim
    real, intent(in)     :: sk, g
    real, intent(out)    :: pos(nx,3), vel(nx,3), acc(nx,3), mas(nx), den(nx), sln(nx), prs(nx), uie(nx), cf(nx), kcf(nx)
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
      spatVarBrdrs21 = 0
      spatVarBrdrs22 = 0.5
    end if
    if (dim.eq.3) then
      spatVarBrdrs31 = 0.
      spatVarBrdrs32 = 0.05
    end if

    parSpacing1 = .01
    parSpacing2 = .01

    shockDensity1 = 1000.
    shockDensity2 = 1000.

    n = 1
    nbnewX = 1
    nbnewY = 1
    nbnewZ = 1

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
            prs(n) = 0
            uie(n) = 0.
            cf(n)  = 0.
            kcf(n) = 10.
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
    ns = n

    allocate(borderX(nbnewX))
    borderX = brdarrX(1:nbnewX)
    allocate(borderY(nbnewY))
    borderY = brdarrY(1:nbnewY)
    allocate(borderZ(nbnewZ))
    borderZ = brdarrZ(1:nbnewZ)
    print *, '#    placed:', n
    print *, '#  border-x:', nbnewX
    print *, '#  border-y:', nbnewY
    print *, '#  border-z:', nbnewZ
  end subroutine tempr_homog01

  subroutine tempr_set_fixed1(A)
    real, intent(out) :: A(ns)
    integer           :: dim

    call get_dim(dim)

    A(borderX) = 0.
    if(dim.gt.1) then
      A(borderY) = 0.
      if(dim.eq.3) then
        A(borderZ) = 0.
      end if
    end if
  end subroutine tempr_set_fixed1

  subroutine tempr_set_fixed3(A)
    real, intent(out) :: A(ns,3)
    integer           :: dim

    call get_dim(dim)

    A(borderX,1) = 0.
    if(dim.gt.1) then
      A(borderY,2) = 0.
      if(dim.eq.3) then
        A(borderZ,3) = 0.
      end if
    end if
  end subroutine tempr_set_fixed3
end module tempr_setup
