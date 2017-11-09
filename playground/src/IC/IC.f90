module IC
  use omp_lib

  use timing,       only: addTime
  use ArrayResize,  only: resize
  use const
  use kernel,       only: get_krad
  use state,        only: get_tasktype,&
                          getkerntype,&
                          getdim,&
                          ginitvar
  use BC
  use initpositions,  only: uniform,&
                            place_close_packed_fcc

  implicit none

  public :: setupIC

  private
  integer(8) :: start=0, finish=0
contains

  subroutine setupIC(n, sk, g, pspc1, pspc2, &
    x, v, dv, mas, den, sln, prs, iu, du, cf, kcf, dcf, ptype)
    integer, allocatable, intent(inout) :: ptype(:)
    real, allocatable, intent(inout), dimension(:,:)  :: x, v, dv, cf, dcf
    real, allocatable, intent(inout), dimension(:,:,:):: kcf
    real, allocatable, intent(inout), dimension(:)    :: mas, den, sln, prs, iu, du
    real, intent(in)     :: sk
    real, intent(inout)  :: pspc1, pspc2, g
    integer, intent(out) :: n

    real                 :: kr, prs1, prs2, rho1, rho2, kcf1, kcf2, cf1, cf2, sp, v0, &
                            brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, period, eA, &
                            cca(3)
    integer              :: i, nb, tt, kt, dim, nptcs, ivt

    call system_clock(start)

    call getkerntype(kt)
    call get_tasktype(tt)
    call get_krad(kr)
    call getdim(dim)
    call ginitvar(ivt)

    if ( kt == 3 ) then
      kr = kr * 2
    end if

    select case (tt)
    case (1)
      ! hydroshock
      ! pspc1 = 0.001
      pspc2 = pspc1*2.
      if (dim == 1) then
        pspc2 = pspc1 * 8.
      end if
      nb = int(kr * sk * 3)
      brdx1 = -.5
      brdx2 = .5
      if ( dim > 1) then
        brdy1 = -pspc2*nb*2
        brdy2 = pspc2*nb*2
      else
        brdy1 = 0.
        brdy2 = 0.
      end if
      if ( dim == 3 ) then
        brdz1 = -pspc2*nb
        brdz2 = pspc2*nb
      else
        brdz1 = 0.
        brdz2 = 0.
      end if
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, x, ptype)
    case (2)
      ! infslb
      nb = 1
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, x, ptype)
    case (3, 9)
      ! heatconduction ! soundwave
      brdx1 = -1.
      brdx2 =  1.
      nptcs = int((brdx2-brdx1)/pspc1)
      pspc1 = merge(0.,(brdx2-brdx1)/nptcs, nptcs == 0)
      pspc1 = pspc2
      nb = int(kr * sk)*2
      if (dim > 1) then
        if (ivt == 4) then
          brdy1 = -1.
          brdy2 =  1.
        else
          brdy1 = -pspc1*nb*2
          brdy2 =  pspc1*nb*2
        end if
      else
        brdy1 = 0.
        brdy2 = 0.
      end if
      if (dim == 3) then
        if (ivt == 4) then
          brdz1 = -1.
          brdz2 =  1.
        else
          brdz1 = -pspc1*nb*2
          brdz2 =  pspc1*nb*2
        end if
      else
        brdz1 = 0.
        brdz2 = 0.
      end if
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, x, ptype)
    case (4)
      ! pheva
      nb = int(kr * sk) + 1
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, x, ptype)
    case(5, 6, 7, 8, 10)
      ! diff-laplace ! diff-graddiv ! diff-artvisc
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
      nb = int(kr * sk)
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, x, ptype)
      ! call place_close_packed_fcc(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, nb, x)
    case default
      print *, 'Task type was not defined in IC.f90: line 125'
      stop
    end select

    n = size(x,dim=2)
    if (n < 2) then
      error stop 'There is only ' // char(n+48) // ' particles in the domain'
    end if

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
    allocate(cf(3,n))
    cf(:,:) = 0.
    allocate(dcf(3,n))
    dcf(:,:) = 0.

    select case (tt)
    case (2, 3)
      allocate(kcf(3,3,n))
      kcf(:,:,:) = 0.
    end select


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
      prs2 = prs1 / 10.
      rho1 = 1.
      rho2 = rho1 / 8.
    case (2)
      ! infslb
      g = 5./3.
      kcf1 = 1.
      kcf2 = 10.
      cf1 = 0.
      cf2 = 1.
    case (3)
      ! heatconduction
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
    case(5, 6, 7, 8, 10)
      ! diff-laplace ! diff-graddiv ! chi-laplace ! chi-graddiv ! diff-artvisc
      rho1 = 1.
    case(9)
      ! soundwave
      rho1 = 1.
      g = 5./3.
      eA = 0.005
    case default
      print *, 'Task type was not defined in IC.f90: line 200'
      stop
    end select

    !-------------------------------!
    !        particles values       !
    !-------------------------------!

    !$omp parallel do default(none)&
    !$omp private(i, sp, cca)&
    !$omp shared(n, pspc1, pspc2, x, sln, den, prs, mas, iu, g, rho1, rho2)&
    !$omp shared(cf, kcf, dim, sk, tt, prs1, prs2, cf1, cf2, kcf1, kcf2)&
    !$omp shared(brdx2, brdx1, brdy2, brdy1, brdz2, brdz1, v0, v, period, ptype)&
    !$omp shared(ivt, eA)
    do i=1,n
      sp = merge(pspc1, pspc2, x(1,i) < 0)
      sln(i) = sk * sp

      select case (tt)
      case (1)
        ! hydroshock
        if (x(1,i) <= 0.) then
          sln(i) = sk * pspc1
          den(i) = rho1
          prs(i) = prs1
          mas(i) = (pspc1**dim) * rho1
          iu(i)  = prs1/(g-1)/rho1
        else
          sln(i) = sk * pspc2
          den(i) = rho2
          prs(i) = prs2
          mas(i) = (pspc2**dim) * rho2
          iu(i)  = prs2/(g-1)/rho2
        end if
      case (2)
        ! infslb
        if (x(1,i) < 0) then
          cf(:,i)  = cf1
          kcf(1,1,i) = kcf1
          kcf(2,2,i) = kcf1
          kcf(3,3,i) = kcf1
          mas(i) = (sp**dim) * rho1
        else
          cf(:,i)  = cf2
          kcf(1,1,i) = kcf2
          kcf(2,2,i) = kcf2
          kcf(3,3,i) = kcf2
          mas(i) = (sp**dim) * rho2
        end if
      case (3, 4)
        mas(i) = (sp**dim) * rho1
        den(i) = rho1
      case(5, 6, 7, 8, 10)
        ! diff-graddiv ! diff-laplace ! chi-laplace
        ! chi-graddiv  ! soundwave ! diff-artvisc
        g = 5./3.
        den(i) = rho1
        mas(i) = (sp**dim) * rho1
        ! v(:,i)  = sin(period*x(:,i))*cos(period*x(:,i)**2)
        v(:,i) = 0.
        if (dim == 1) then
          ! x
          ! v(1,i) = x(1,i)
          ! x sinx
          v(1,i) = sin(x(1,i)) * x(1,i)
          ! sin
          ! v(1,i) = sin(x(1,i))
        elseif ( dim == 2 ) then
          ! v(1,i) = x(1,i)*x(2,i)
          ! v(2,i) = x(1,i)*x(2,i)
          ! (y six, x^2 siny)
          v(1,i) = sin(x(1,i)) * x(2,i)
          v(2,i) = sin(x(2,i)) * x(1,i) * x(1,i)
          ! sin
          ! v(1,i) = sin(x(1,i))
          ! v(2,i) = sin(x(2,i))
        elseif ( dim == 3 ) then
          ! v(1,i) = x(1,i)*x(2,i)*x(3,i)
          ! v(2,i) = x(1,i)*x(2,i)*x(3,i)
          ! v(3,i) = x(1,i)*x(2,i)*x(3,i)
          ! (y six, z^2 siny, x^3 sinz)
          v(1,i) = sin(x(1,i)) * x(2,i)
          v(2,i) = sin(x(2,i)) * x(3,i) * x(3,i)
          v(3,i) = sin(x(3,i)) * x(1,i) * x(1,i) * x(1,i)
          ! sin
          ! v(1,i) = sin(x(1,i))
          ! v(2,i) = sin(x(2,i))
          ! v(3,i) = sin(x(3,i))
        end if
      case(9)
        den(i) = rho1 * (1. + eA * sin(pi * x(1,i)))
        mas(i) = (sp**dim) * den(i)
        v(:,i) = 0.
        v(1,i) = eA*sin(pi*(x(1,i)))
      case default
        print *, 'Task type was not defined in IC.f90: line 300'
        stop
      end select

      select case (ivt)
      case(-1)
      case(1)
        ! isotropic-sinxsinysinz
        kcf(1,1,i) = kcf1
        kcf(2,2,i) = kcf1
        kcf(3,3,i) = kcf1
        ! if (x(1,i) < 0) then
        !   kcf(i) = 1
        ! else
        !   kcf(i) = 10
        ! end if
        cf(:, i) = 0
        if (ptype(i) /= 0) then
          if ( dim == 1) then
            cf(1, i)  = sin(pi * (x(1,i) - brdx1) / abs(brdx2-brdx1))
          elseif ( dim == 2 ) then
            cf(1, i)  = sin(pi * (x(1,i) - brdx1) / abs(brdx2-brdx1)) * &
                     sin(pi * (x(2,i) - brdy1) / abs(brdy2-brdy1))
          elseif ( dim == 3 ) then
            cf(1, i)  = sin(pi * (x(1,i) - brdx1) / abs(brdx2-brdx1)) * &
                     sin(pi * (x(2,i) - brdy1) / abs(brdy2-brdy1)) * &
                     sin(pi * (x(3,i) - brdz1) / abs(brdz2-brdz1))
          end if
        end if
      case(2)
        ! orthotropic-sinxsinysinz
        kcf(1,1,i) = 10
        kcf(2,2,i) = 1
        kcf(3,3,i) = 0.1
        cf(:, i) = 0
        if (ptype(i) /= 0) then
          if ( dim == 1) then
            cf(1, i)  = sin(pi * (x(1,i) - brdx1) / abs(brdx2-brdx1))
          elseif ( dim == 2 ) then
            cf(1, i)  = sin(pi * (x(1,i) - brdx1) / abs(brdx2-brdx1)) * &
                     sin(pi * (x(2,i) - brdy1) / abs(brdy2-brdy1))
          elseif ( dim == 3 ) then
            cf(1, i)  = sin(pi * (x(1,i) - brdx1) / abs(brdx2-brdx1)) * &
                     sin(pi * (x(2,i) - brdy1) / abs(brdy2-brdy1)) * &
                     sin(pi * (x(3,i) - brdz1) / abs(brdz2-brdz1))
          end if
        end if
      case(3)
        ! shock12
        kcf(1,1,i) = 1.
        kcf(1,2,i) = 0.
        kcf(1,3,i) = 0.
        kcf(2,1,i) = 0.
        kcf(2,2,i) = 0.
        kcf(2,3,i) = 0.
        kcf(3,1,i) = 0.
        kcf(3,2,i) = 0.
        kcf(3,3,i) = 0.
        ! cf(:, i) = sin(x(:,i))
        if (x(1,i) < 0.) then
          cf(1,i) = 1
        else if (x(1,i) > 0.) then
          cf(1,i) = 2
        else
          cf(1,i) = 1.5
        end if
      case(4)
        ! gaussian-pulse
        kcf(1,1,i) = 1.
        kcf(1,2,i) = 0.
        kcf(1,3,i) = 0.
        kcf(2,1,i) = 0.
        kcf(2,2,i) = 0.
        kcf(2,3,i) = 0.
        kcf(3,1,i) = 0.
        kcf(3,2,i) = 0.
        kcf(3,3,i) = 0.
        ! cf(:, i) = sin(x(:,i))
        cf(1,i) = (2.*pi)**(-dim/2.)/(0.05**2)**(dim/2.)*&
                  exp(-0.5*(x(1,i)*x(1,i) + x(2,i)*x(2,i) + x(3,i)*x(3,i))/(0.05**2))
      case(5)
        ! gaussian-ring
        kcf(1,1,i) = 0.
        kcf(1,2,i) = 0.
        kcf(1,3,i) = 0.
        kcf(2,1,i) = 0.
        kcf(2,2,i) = 1.
        kcf(2,3,i) = 0.
        kcf(3,1,i) = 0.
        kcf(3,2,i) = 0.
        kcf(3,3,i) = 0.

        cca(1) = sqrt(x(1,i)*x(1,i) + x(2,i)*x(2,i))
        cca(2) = atan(x(2,i),x(1,i))
        cca(3) = x(3,i)

        ! exp[−(1/2)[(r−r0)^2/δr0^2 + φ^2/δφ0^2]]
        ! δr0 = 0.05 and r0 = 0.3 define a Gaussian ring
        !  at radius r0 of width δr, and δφ0 = 0.5
        ! is an initial Gaussian spread about φ=0intheφˆdirection
        cf(1,i) = exp(-0.5*((cca(1) - 0.3)**2/(0.05*0.05) + cca(2)*cca(2)/(0.5*0.5)))
      end select
    end do
    !$omp end parallel do
    call system_clock(finish)
    call addTime(' ic', finish - start)
  end subroutine
end module
