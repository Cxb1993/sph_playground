module IC
  use utils
  use kernel
  use BC
  use uniform
  use semiuniform

  implicit none

  public :: setup

  private
  real, parameter     :: pi = 4.*atan(1.)

contains

  subroutine setup(tt, kt, dim, n, sk, g, cv, pspc1, pspc2, &
    pos, vel, acc, mas, den, sln, prs, iu, du, cf, kcf, dcf)
    real, allocatable, intent(inout) :: pos(:,:), vel(:,:), acc(:,:), mas(:), den(:), sln(:), &
                                        prs(:), iu(:), du(:), cf(:), kcf(:), dcf(:)
    character(len=*), intent(in) :: tt, kt
    integer, intent(in)  :: dim
    real, intent(in)     :: sk, cv
    real, intent(inout)  :: pspc1, pspc2, g
    integer, intent(out) :: n
    real                 :: kr, prs1, prs2, rho1, rho2, kcf1, kcf2, cf1, cf2, sp, v0, &
                            brdx1, brdx2, brdy1, brdy2, brdz1, brdz2
    integer              :: i, nb

    call set_dim(dim)
    call set_tasktype(tt)
    call set_kerntype(kt)
    call get_krad(kr)

    nb = 0
    brdx1 = -1.
    brdx2 = 1.
    if (dim > 1) then
      brdy1 = - int(kr * sk) * pspc1 * 2
      brdy2 =   int(kr * sk) * pspc1 * 2
    else
      brdy1 = 0.
      brdy2 = 0.
    end if
    if (dim == 3) then
      brdz1 = - int(kr * sk) * pspc1 * 2
      brdz2 =   int(kr * sk) * pspc1 * 2
    else
      brdz1 = 0.
      brdz2 = 0.
    end if

    select case (tt)
    case ('hydroshock')
      brdx1 = -.5
      brdx2 = .5
      nb = 4
      pspc1 = 0.001
      pspc2 = 0.008
      call make_semiuniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos)
    case ('infslb')
      nb = 1
      call make_uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos)
    case ('hc-sinx')
      nb = int(kr * sk) + 1
      call make_uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos)
    case ('pheva')
      nb = int(kr * sk) + 1
      call make_uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos)
    case default
      print *, 'Task type was not defined in IC border stage'
      stop
    end select

    n = size(pos,dim=2)

    allocate(vel(3,n))
    vel(:,:) = 0.
    allocate(acc(3,n))
    acc(:,:) = 0.
    allocate(mas(n))
    mas(:) = 0.
    allocate(den(n))
    den(:) = 0.
    allocate(sln(n))
    sln(:) = 0.
    allocate(prs(n))
    prs(:) = 0.
    allocate(iu(n))
    iu(:) = 0.
    allocate(du(n))
    du(:) = 0.
    allocate(cf(n))
    cf(:) = 0.
    allocate(kcf(n))
    kcf(:) = 0.
    allocate(dcf(n))
    dcf(:) = 0.
    !
    ! common values
    !
    rho1 = 1.
    rho2 = 1.
    prs1 = 1.
    prs2 = 1.
    kcf1 = 1.
    kcf2 = 1.
    cf1  = 0.

    select case (tt)
    case ('hydroshock')
      g = 1.4
      prs1 = 1.
      prs2 = 0.1
      rho1 = 1.
      rho2 = 0.125
    case ('infslb')
      g = 1.4
      kcf1 = 1.
      kcf2 = 10.
      cf1 = 0.
      cf2 = 1.
    case ('hc-sinx')
    case ('pheva')
      g    = 5/3.
      rho1 = 2.
      cf1  = 0.5
      kcf1 = 1e-1
      v0   = 1e-4
      prs1 = 1.
    case default
      print *, 'Task type was not defined in IC state stage'
      stop
    end select

    do i=1,n
      sp = merge(pspc1, pspc2, pos(1,i) < 0)
      sln(i) = sk * sp

      select case (tt)
      case ('hydroshock')
        if (pos(1,i) < 0) then
          den(i) = rho1
          prs(i) = prs1
          mas(i) = (sp**dim) * rho1
        else
          den(i) = rho2
          prs(i) = prs2
          mas(i) = (sp**dim) * rho2
        end if
        iu(i) = merge(prs1/(g-1)/rho1, prs2/(g-1)/rho2, pos(1,i) < 0)
      case ('infslb')
        if (pos(1,i) < 0) then
          cf(i)  = cf1
          kcf(i) = kcf1
          mas(i) = (sp**dim) * rho1
        else
          cf(i)  = cf2
          kcf(i) = kcf2
          mas(i) = (sp**dim) * rho2
        end if
        iu(i) = cf(i) / cv
      case ('hc-sinx')
        mas(i) = (sp**dim) * rho1
        den(i) = rho1
        kcf(i) = kcf1
        prs(i) = prs1
        cf(i)  = sin(2 * pi * (pos(1,i) - brdx2) / abs(brdx2-brdx1))
        iu(i) = cf(i) / cv
      case ('pheva')
        cf(i)  = cf1
        den(i) = rho1 * (1 + v0 * sin(2 * pi * (pos(1,i) - brdx2) / abs(brdx2-brdx1)))
        mas(i) = (sp**dim) * den(i)
        vel(1,i) = v0 * sin(2 * pi * (pos(1,i) - brdx2) / abs(brdx2-brdx1))
        kcf(i) = kcf1
        prs(i) = prs1
        iu(i)  = prs1/(g-1)/(1-cf(i))/rho1
      case default
        print *, 'Task type was not defined in IC state stage'
        stop
      end select
    end do
  end subroutine setup
end module IC
