module IC
  use timing, only: addTime
  use utils
  use const
  use kernel, only: get_tasktype,&
                    get_kerntype,&
                    get_krad,&
                    get_dim
  use BC
  use initpositions,  only: uniform,&
                            place_close_packed_fcc
  use semiuniform

  implicit none

  public :: setupIC

  private
  integer(8) :: start=0, finish=0
contains

  subroutine setupIC(n, sk, g, cv, pspc1, pspc2, &
    x, v, dv, mas, den, sln, prs, iu, du, cf, kcf, dcf, ptype)
    integer, allocatable, intent(inout) :: ptype(:)
    real, allocatable, intent(inout) :: x(:,:), v(:,:), dv(:,:), mas(:), den(:), sln(:), &
                                        prs(:), iu(:), du(:), cf(:), kcf(:), dcf(:)
    real, intent(in)     :: sk, cv
    real, intent(inout)  :: pspc1, pspc2, g
    integer, intent(out) :: n

    real                 :: kr, prs1, prs2, rho1, rho2, kcf1, kcf2, cf1, cf2, sp, v0, &
                            brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, period
    integer              :: i, nb, tt, kt, dim, nptcs

    call system_clock(start)

    call get_kerntype(kt)
    call get_tasktype(tt)
    call get_krad(kr)
    call get_dim(dim)

    if ( kt == 3 ) then
      kr = kr * 2
    end if

    select case (tt)
    case (1)
      ! hydroshock
      brdx1 = -.5
      brdx2 = .5
      nb = 4
      pspc1 = 0.001
      pspc2 = 0.008
      call make_semiuniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, x)
    case (2)
      ! infslb
      nb = 1
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, x, ptype)
    case (3)
      ! hc-sinx ! chi-laplace
      brdx1 = -1.
      brdx2 =  1.
      nptcs = int((brdx2-brdx1)/pspc1)
      pspc1 = merge(0.,(brdx2-brdx1)/nptcs, nptcs == 0)
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
      nb = int(kr * sk) + 1
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, x, ptype)
    case (4)
      ! pheva
      nb = int(kr * sk) + 1
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, x, ptype)
    case(5, 6, 7, 8)
      ! diff-laplace ! diff-graddiv
      period = pi
      brdx1 = -1.*period
      brdx2 = 1.*period
      if (dim > 1) then
        brdy1 = -1.*period
        brdy2 =  1.*period
      else
        brdy1 = 0.
        brdy2 = 0.
      end if
      if (dim == 3) then
        brdz1 = -1.*period
        brdz2 =  1.*period
      else
        brdz1 = 0.
        brdz2 = 0.
      end if
      nb = int(kr * sk) + 1
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, x, ptype)
      ! call place_close_packed_fcc(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, nb, x)
    case default
      print *, 'Task type was not defined in IC border stage'
      stop
    end select

    ! print *, brdy1, brdy2, pspc1, x(:,size(x,2))
    ! read *
    !
    n = size(x,dim=2)

    allocate(v(3,n))
    v(:,:) = 0.
    allocate(dv(3,n))
    dv(:,:) = 0.
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

    !--------------------
    ! common values
    !--------------------
    rho1 = 0.
    rho2 = 0.
    prs1 = 0.
    prs2 = 0.
    kcf1 = 0.
    kcf2 = 0.
    cf1  = 0.

    select case (tt)
    case (1)
      ! hydroshock
      g = 1.4
      prs1 = 1.
      prs2 = 0.1
      rho1 = 1.
      rho2 = 0.125
    case (2)
      ! infslb
      g = 1.4
      kcf1 = 1.
      kcf2 = 10.
      cf1 = 0.
      cf2 = 1.
    case (3)
      ! hc-sinx ! chi-laplace
      rho1 = 1.
      kcf1 = 1.
    case (4)
      ! pheva
      g    = 5/3.
      rho1 = 2.
      cf1  = 0.5
      kcf1 = 1e-1
      v0   = 1e-4
      prs1 = 1.
    case(5, 6, 7, 8)
      ! diff-laplace ! diff-graddiv
      rho1 = 1.
      period = 1.
    case default
      print *, 'Task type was not defined in IC state stage'
      stop
    end select

    !-----------------
    ! particles values
    !-----------------
    do i=1,n
      sp = merge(pspc1, pspc2, x(1,i) < 0)
      sln(i) = sk * sp

      select case (tt)
      case (1)
        ! hydroshock
        if (x(1,i) < 0) then
          den(i) = rho1
          prs(i) = prs1
          mas(i) = (sp**dim) * rho1
        else
          den(i) = rho2
          prs(i) = prs2
          mas(i) = (sp**dim) * rho2
        end if
        iu(i) = merge(prs1/(g-1)/rho1, prs2/(g-1)/rho2, x(1,i) < 0)
      case (2)
        ! infslb
        if (x(1,i) < 0) then
          cf(i)  = cf1
          kcf(i) = kcf1
          mas(i) = (sp**dim) * rho1
        else
          cf(i)  = cf2
          kcf(i) = kcf2
          mas(i) = (sp**dim) * rho2
        end if
        iu(i) = cf(i) / cv
      case (3)
        ! hc-sinx
        mas(i) = (sp**dim) * rho1
        den(i) = rho1
        kcf(i) = kcf1
        prs(i) = prs1
        cf(i)  = sin(2 * pi * (x(1,i) - brdx2) / abs(brdx2-brdx1))
        iu(i) = cf(i) / cv
      case (4)
        ! pheva
        cf(i)  = cf1
        den(i) = rho1 * (1 + v0 * sin(2 * pi * (x(1,i) - brdx2) / abs(brdx2-brdx1)))
        mas(i) = (sp**dim) * den(i)
        v(1,i) = v0 * sin(2 * pi * (x(1,i) - brdx2) / abs(brdx2-brdx1))
        kcf(i) = kcf1
        prs(i) = prs1
        iu(i)  = prs1/(g-1)/(1-cf(i))/rho1
      case(5, 6, 7, 8)
        ! diff-graddiv ! diff-laplace ! chi-laplace ! chi-graddiv
        den(i) = rho1
        mas(i) = (sp**dim) * rho1
        ! v(:,i)  = sin(period*x(:,i))*cos(period*x(:,i)**2)
        v(:,i) = 0.
        if (dim == 1) then
          ! x sinx
          ! v(1,i) = sin(period*x(1,i)) * x(1,i)
          ! sin
          v(1,i) = sin(period*x(1,i))
          ! v(1,i) = x(1,i)
        elseif ( dim == 2 ) then
          ! v(1,i) = x(1,i)*x(2,i)
          ! v(2,i) = x(1,i)*x(2,i)
          ! ysix xsiny
          ! v(1,i) = sin(period*x(1,i))* x(2,i)
          ! v(2,i) = sin(period*x(2,i))* x(1,i)
          ! sin
          v(1,i) = sin(period*x(1,i))
          v(2,i) = sin(period*x(2,i))
        elseif ( dim == 3 ) then
          ! v(1,i) = x(1,i)*x(2,i)*x(3,i)
          ! v(2,i) = x(1,i)*x(2,i)*x(3,i)
          ! v(3,i) = x(1,i)*x(2,i)*x(3,i)
          ! ysix zsiny xsinz
          ! v(1,i) = sin(period*x(1,i))* x(2,i)
          ! v(2,i) = sin(period*x(2,i))* x(3,i)
          ! v(3,i) = sin(period*x(3,i))* x(1,i)
          ! sin
          v(1,i) = sin(period*x(1,i))
          v(2,i) = sin(period*x(2,i))
          v(3,i) = sin(period*x(3,i))
        end if
      case default
        print *, 'Task type was not defined in IC state stage'
        stop
      end select
    end do
    call system_clock(finish)
    call addTime(' ic', finish - start)
  end subroutine setupIC
end module IC
