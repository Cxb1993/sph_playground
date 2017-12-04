module IC
  use omp_lib

  use timing,       only: addTime
  use ArrayResize,  only: resize
  use const
  use kernel,       only: get_krad, &
                          get_w
  use state,        only: get_tasktype,&
                          getkerntype,&
                          getdim,&
                          ginitvar,&
                          gcoordsys,&
                          sorigin
  use BC
  use initpositions,  only: uniform,&
                            place_close_packed_fcc

  use neighboursearch,  only: getneighbours, &
                              getNeibListL1, &
                              getNeibListL2, &
                              findneighboursKDT
  implicit none

  public :: setupIC

  private
  integer(8) :: start=0, finish=0
contains

  subroutine setupIC(n, sk, g, pspc1, pspc2, &
    pos, v, dv, mas, den, sln, prs, iu, du, cf, kcf, dcf, ptype)
    integer, allocatable, intent(inout) :: ptype(:)
    real, allocatable, intent(inout), dimension(:,:)  :: pos, v, dv, cf, dcf
    real, allocatable, intent(inout), dimension(:,:,:):: kcf
    real, allocatable, intent(inout), dimension(:)    :: mas, den, sln, prs, iu, du
    real, intent(in)     :: sk
    real, intent(inout)  :: pspc1, pspc2, g
    integer, intent(out) :: n

    real                 :: kr, prs1, prs2, rho1, rho2, kcf1, kcf2, cf1, cf2, sp, v0, &
                            brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, period, eA, &
                            cca(3), qmatr(3,3), qtmatr(3,3), w
    integer              :: i, lj, j, nb, tt, kt, dim, nptcs, ivt, cs
    integer, allocatable :: nlista(:)
    integer(8)           :: t0

    call system_clock(start)

    call getkerntype(kt)
    call get_tasktype(tt)
    call get_krad(kr)
    call getdim(dim)
    call ginitvar(ivt)
    call gcoordsys(cs)

    if ( kt == 3 ) then
      kr = kr * 2
    end if

    !
    !
    ! need to get rid of it
    !
    !
    select case (tt)
    case (1, 2, 9, 4)
      ! mooved to ivt check below
    case (3)
      ! heatconduction
      brdx1 = -1.
      brdx2 =  1.
      nptcs = int((brdx2-brdx1)/pspc1)
      pspc1 = merge(0.,(brdx2-brdx1)/nptcs, nptcs == 0)
      pspc2 = pspc1
      nb = int(kr * sk)*2
      if (dim > 1) then
        if ((ivt == 4).or.(ivt == 5)) then
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
        if ((ivt == 4).or.(ivt == 5)) then
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
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos, ptype, randomise=0.)
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
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos, ptype, randomise=0.)
      ! call place_close_packed_fcc(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, nb, pos)
    case default
      print *, 'Task type was not defined in IC.f90: line 170'
      stop
    end select
    !
    !
    ! need to get rid of it
    !
    !

    select case(ivt)
    case(-1, 1, 2, 3, 4, 5)
      print*, "Not set yet. FIX ME."
    case (6)
      ! soundwave
      rho1 = 1.
      g = 5./3.
      eA = 0.005
      prs1 = 1.

      brdx1 = -1.
      brdx2 =  1.
      nptcs = int((brdx2-brdx1)/pspc1)
      pspc1 = merge(0.,(brdx2-brdx1)/nptcs, nptcs == 0)
      pspc2 = pspc1
      nb = int(kr * sk * 2.)
      if (dim > 1) then
        brdy1 = -pspc1*nb*2
        brdy2 =  pspc1*nb*2
      else
        brdy1 = 0.
        brdy2 = 0.
      end if
      if (dim == 3) then
        brdz1 = -pspc1*nb*2
        brdz2 =  pspc1*nb*2
      else
        brdz1 = 0.
        brdz2 = 0.
      end if
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos, ptype, randomise=0.)
    case (7)
      ! hydroshock
      g = 1.4
      prs1 = 1.
      prs2 = prs1 / 10.
      rho1 = 1.
      rho2 = rho1 / 8.
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
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos, ptype, randomise=0.)
    case (8)
      ! alfvenwave
      rho1 = 1.
      g    = 5./3.
      eA   = 0.1
      prs1 = 0.1

      brdx1 = -1.5
      brdx2 =  1.5
      nptcs = int((brdx2-brdx1)/pspc1)
      pspc1 = merge(0.,(brdx2-brdx1)/nptcs, nptcs == 0)
      pspc2 = pspc1
      nb = int(kr * sk*1.5)
      if (dim > 1) then
        brdy1 = -0.5
        brdy2 =  0.5
      else
        brdy1 = 0.
        brdy2 = 0.
      end if
      if (dim == 3) then
        brdz1 = -0.5
        brdz2 =  0.5
      else
        brdz1 = 0.
        brdz2 = 0.
      end if
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos, ptype, randomise=0.)
    case default
      print *, 'Problem was not set in IC.f90: line 215'
      stop
    end select

    n = size(pos,dim=2)
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
    allocate(kcf(3,3,n))
    kcf(:,:,:) = 0.

    select case (tt)
    case (1, 2, 4, 9)
      ! mover to ivt check above
    case (3)
      ! heatconduction
      rho1 = 1.
      kcf1 = 1.
    case(5, 6, 7, 8, 10)
      ! diff-laplace ! diff-graddiv ! chi-laplace ! chi-graddiv ! diff-artvisc
      g    = 5/3.
      rho1 = 1.
    case default
      print *, 'Task type was not defined in IC.f90: line 200'
      stop
    end select

    !-------------------------------!
    !        particles values       !
    !-------------------------------!

    ! sln(:) = sk * pspc1
    ! call findneighboursKDT(ptype, pos, sln)

    !$omp parallel do default(none)&
    !$omp private(i, sp, cca, qmatr, qtmatr, t0, lj, j, w, nlista)&
    !$omp shared(n, pspc1, pspc2, pos, sln, den, prs, mas, iu, g, rho1, rho2)&
    !$omp shared(cf, kcf, dim, sk, tt, prs1, prs2, cf1, cf2, kcf1, kcf2)&
    !$omp shared(brdx2, brdx1, brdy2, brdy1, brdz2, brdz1, v0, v, period, ptype)&
    !$omp shared(ivt, eA, kt, cs)
    do i=1,n
      sp = merge(pspc1, pspc2, pos(1,i) < 0)
      select case (tt)
      case (1, 2)
        ! mooved to ivt check below
      case (3, 4)
        sln(i) = sk * sp
        mas(i) = (sp**dim) * rho1
        den(i) = rho1
        ! call getneighbours(i, pos, sln, nlista, t0)
        ! mas(i) = 0.
        ! do lj = 1, size(nlista)
        !   j = nlista(lj)
        !   call get_w(sqrt(dot_product(pos(:,i) - pos(:,j),pos(:,i) - pos(:,j))), sln(i), w)
        !   mas(i) = mas(i) + w / rho1
        ! end do
        ! call get_w(0., sln(i), w)
        ! mas(i) = mas(i) + w / rho1
        ! mas(i) = 1./mas(i)
        ! den(i) = mas(i)/(sp**dim)
        ! print*, '# ', i, '-------------------------'
        ! print*, 'm_i_numer = ', mas(i)
        ! print*, 'm_i_theor = ', (sp**dim) * rho1
        ! print*, 'accuracy  % ', ((sp**dim) * rho1)/mas(i)
        ! print*, ''
        ! print*, 'r_i_numer = ', den(i)
        ! print*, 'r_i_theor = ', rho1
        ! print*, 'accuracy  % ', rho1/den(i)
        ! print*, ''
        ! print*, 'n_i_numer = ', size(nlista)
        ! print*, 'n_i_theor = ', pi*(2.*sk)**2
        ! print*, 'accuracy  % ', size(nlista)/(pi*(2.*sk)**2)
        ! read*
        ! den(i) = rho1
      case(5, 6, 7, 8, 10)
        ! diff-graddiv ! diff-laplace ! chi-laplace
        ! chi-graddiv  ! soundwave ! diff-artvisc
        sln(i) = sk * sp
        den(i) = rho1
        mas(i) = (sp**dim) * rho1
        ! v(:,i)  = sin(period*pos(:,i))*cos(period*pos(:,i)**2)
        v(:,i) = 0.
        if (dim == 1) then
          ! pos
          ! v(1,i) = pos(1,i)
          ! pos sinpos
          v(1,i) = sin(pos(1,i)) * pos(1,i)
          ! sin
          ! v(1,i) = sin(pos(1,i))
        elseif ( dim == 2 ) then
          ! v(1,i) = pos(1,i)*pos(2,i)
          ! v(2,i) = pos(1,i)*pos(2,i)
          ! (y six, x^2 siny)
          v(1,i) = sin(pos(1,i)) * pos(2,i)
          v(2,i) = sin(pos(2,i)) * pos(1,i) * pos(1,i)
          ! sin
          ! v(1,i) = sin(pos(1,i))
          ! v(2,i) = sin(pos(2,i))
        elseif ( dim == 3 ) then
          ! v(1,i) = pos(1,i)*pos(2,i)*pos(3,i)
          ! v(2,i) = pos(1,i)*pos(2,i)*pos(3,i)
          ! v(3,i) = pos(1,i)*pos(2,i)*pos(3,i)
          ! (y six, z^2 siny, x^3 sinz)
          v(1,i) = sin(pos(1,i)) * pos(2,i)
          v(2,i) = sin(pos(2,i)) * pos(3,i) * pos(3,i)
          v(3,i) = sin(pos(3,i)) * pos(1,i) * pos(1,i) * pos(1,i)
          ! sin
          ! v(1,i) = sin(pos(1,i))
          ! v(2,i) = sin(pos(2,i))
          ! v(3,i) = sin(pos(3,i))
        end if
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
        ! if (pos(1,i) < 0) then
        !   kcf(i) = 1
        ! else
        !   kcf(i) = 10
        ! end if
        cf(:, i) = 0
        if (ptype(i) /= 0) then
          if ( dim == 1) then
            cf(1, i)  = sin(pi * (pos(1,i) - brdx1) / abs(brdx2-brdx1))
          elseif ( dim == 2 ) then
            cf(1, i)  = sin(pi * (pos(1,i) - brdx1) / abs(brdx2-brdx1)) * &
                     sin(pi * (pos(2,i) - brdy1) / abs(brdy2-brdy1))
          elseif ( dim == 3 ) then
            cf(1, i)  = sin(pi * (pos(1,i) - brdx1) / abs(brdx2-brdx1)) * &
                     sin(pi * (pos(2,i) - brdy1) / abs(brdy2-brdy1)) * &
                     sin(pi * (pos(3,i) - brdz1) / abs(brdz2-brdz1))
          end if
        end if
      case(3)
        ! shock12
        kcf(1,1,i) = 0.
        kcf(1,2,i) = 0.
        kcf(1,3,i) = 0.
        kcf(2,1,i) = 0.
        kcf(2,2,i) = 1.
        kcf(2,3,i) = 0.
        kcf(3,1,i) = 0.
        kcf(3,2,i) = 0.
        kcf(3,3,i) = 0.
        ! cf(:, i) = sin(pos(:,i))
        if (pos(1,i) < 0.) then
          cf(1,i) = 1.
        else if (pos(1,i) > 0.) then
          cf(1,i) = 2.
        else
          cf(1,i) = 1.5
        end if
      case(4)
        ! pulse
        kcf(1,1,i) = 1.
        kcf(1,2,i) = 0.
        kcf(1,3,i) = 0.
        kcf(2,1,i) = 0.
        kcf(2,2,i) = 0.
        kcf(2,3,i) = 0.
        kcf(3,1,i) = 0.
        kcf(3,2,i) = 0.
        kcf(3,3,i) = 0.
        ! cf(:, i) = sin(pos(:,i))
        cf(1,i) = (2.*pi)**(-dim/2.)/(0.1**2)**(dim/2.)*&
                  exp(-0.5*(pos(1,i)*pos(1,i) + pos(2,i)*pos(2,i) + pos(3,i)*pos(3,i))/(0.1**2))
      case(5)
        ! ring
        ! initially in cylindric CS, but will be transphered in other system few lines lower
        kcf(1,1,i) = 0.
        kcf(1,2,i) = 0.
        kcf(1,3,i) = 0.
        kcf(2,1,i) = 0.
        kcf(2,2,i) = 1.
        kcf(2,3,i) = 0.
        kcf(3,1,i) = 0.
        kcf(3,2,i) = 0.
        kcf(3,3,i) = 0.

        cca(1) = sqrt(pos(1,i)*pos(1,i) + pos(2,i)*pos(2,i))
        cca(2) = atan(pos(2,i),pos(1,i))
        cca(3) = pos(3,i)

        if ( cs == 1 ) then
          ! cartesian
          qmatr(1,1) = cos(cca(2))
          qmatr(1,2) = -sin(cca(2))
          qmatr(1,3) = 0.
          qmatr(2,1) = sin(cca(2))
          qmatr(2,2) = cos(cca(2))
          qmatr(2,3) = 0.
          qmatr(3,1) = 0.
          qmatr(3,2) = 0.
          qmatr(3,3) = 1.

          qtmatr = transpose(qmatr)

          kcf(:,:,i) = matmul(qmatr(:,:), kcf(:,:,i))
          kcf(:,:,i) = matmul(kcf(:,:,i), qtmatr(:,:))
        end if
        ! exp[−(1/2)[(r−r0)^2/δr0^2 + φ^2/δφ0^2]]
        ! δr0 = 0.05 and r0 = 0.3 define a Gaussian ring
        !  at radius r0 of width δr, and δφ0 = 0.5
        ! is an initial Gaussian spread about φ=0intheφˆdirection
        cf(1,i) = exp(-0.5*((cca(1) - 0.3)**2/(0.05*0.05) + cca(2)*cca(2)/(0.5*0.5)))
      case (6)
        ! soundwave
        sln(i) = sk * sp
        den(i) = rho1 * (1. + eA * sin(pi * pos(1,i)))
        mas(i) = (sp**dim) * den(i)
        v(1,i) = eA*sin(pi*(pos(1,i)))
        prs(i) = prs1
        iu(i)  = prs1/(g-1)/rho1
      case (7)
        ! hydroshock
        if (pos(1,i) <= 0.) then
          sln(i) = sk * sp
          den(i) = rho1
          prs(i) = prs1
          mas(i) = (pspc1**dim) * rho1
          iu(i)  = prs1/(g-1)/rho1
        else
          sln(i) = sk * sp
          den(i) = rho2
          prs(i) = prs2
          mas(i) = (pspc2**dim) * rho2
          iu(i)  = prs2/(g-1)/rho2
        end if
      case (8)
        ! alfvenwave
        sln(i) = sk * sp
        mas(i) = (sp**dim) * rho1

        den(i) = rho1
        prs(i) = prs1
        iu(i)  = prs1/(g-1)/rho1

        cca(1) = (pos(1,i) + 2.*pos(2,i) + 2.*pos(3,i))/3.

        v(1,i) = 0.
        v(2,i) = eA*sin(2.*pi*cca(1))
        v(3,i) = eA*cos(2.*pi*cca(1))
        cf(1,i) = 1.
        cf(2,i) = v(2,i)
        cf(3,i) = v(3,i)

        if (dim == 1) then
          v(2:3,i)   = 0.
          cf(2:3,i)  = 0.
        elseif (dim == 2) then
          v(3,i)    = 0.
          cf(3, i)  = 0.
        end if
      end select
    end do
    !$omp end parallel do
    call system_clock(finish)
    call addTime(' ic', finish - start)
  end subroutine
end module
