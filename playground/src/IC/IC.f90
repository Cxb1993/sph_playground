module IC
  use omp_lib

  use errprinter,       only: error
  use timing,           only: addTime
  use ArrayResize,      only: resize
  use const
  use kernel,           only: get_krad, &
                              get_w
  use state,            only: get_equations,&
                              getddwtype,&
                              getdim,&
                              ginitvar,&
                              gcoordsys,&
                              setdiffisotropic, &
                              setdiffconductivity, &
                              setmhdmagneticpressure,&
                              gethfac,&
                              getresolution,&
                              getspacing,&
                              getPartNumber,&
                              setPartNumber,&
                              setGamma
  use BC,               only: createFixedBorders, &
                              getPeriodPartNumber,&
                              createPeriodicBorder,&
                              findInsideBorderParticles,&
                              clearPeriodicParticles
  use placeUniform,     only: uniformV4
  use placeClose,       only: closepacked
  use placeRandom,      only: random
  use neighboursearch,  only: getneighbours, &
                              getNeibListL1, &
                              getNeibListL2, &
                              findneighboursKDT
  use stretchmap,       only: set_density_profile
  use rhofuncs
  use circuit1, only: c1
  use circuit2, only: c2
  implicit none

  public :: setupV2

  private
  integer(8) :: start=0, finish=0
contains

  subroutine setupV2(n, cv, store)
    real, allocatable, intent(inout) :: store(:,:)
    real, intent(in)        :: cv
    integer, intent(inout)  :: n

    real :: &
      ra(3), kr, prs1, prs2, rho1, rho2, sp, &
      brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, eA, &
      cca(3), qmatr(3,3), pspc2, theta, phi, bordersize,&
      hfac, pspc1, gamma, lx, ly, lz
    integer :: &
      i, nb, tt, kt, dim, ivt, cs, d2null, d3null, &
      rpn, fpn, ppn, resol
    ! integer, allocatable :: nlista(:)

    call system_clock(start)
    call getddwtype(kt)
    call get_equations(tt)
    call get_krad(kr)
    call getdim(dim)
    call ginitvar(ivt)
    call gcoordsys(cs)
    call gethfac(hfac)
    call getresolution(resol)
    call getspacing(pspc1)

    d2null = 1
    d3null = 1
    if (dim == 1) then
      d2null = 0
      d3null = 0
    else if (dim == 2) then
      d3null = 0
    end if

    call setmhdmagneticpressure(1.)
    call setdiffisotropic("yes")
    call setdiffconductivity(1.)

    nb = int(kr*hfac)+1

    select case(ivt)
    case(ett_sin3)
      print*, "Not set yet. FIX ME. IC.f90. line 103."
      stop
    case (ett_shock12)
      call setmhdmagneticpressure(1.)
      call setdiffisotropic("yes")
      call setdiffconductivity(1.)
      rho1 = 1.
      prs1 = 1.
      gamma = 5./3.

      brdx1 = -1.
      brdx2 =  1.
      if (resol == 0) then
        resol = int((brdx2-brdx1)/pspc1)
      end if
      pspc1 = (brdx2-brdx1)/resol
      pspc2 = pspc1
      brdy1 = -pspc1*nb*2*d2null
      brdy2 =  pspc1*nb*2*d2null
      brdz1 = -pspc1*nb*2*d3null
      brdz2 =  pspc1*nb*2*d3null
      bordersize = nb*pspc1
      call uniformV4(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, bordersize, pspc1, store, padding=0.5)
      ! call closepacked(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, &
      !   bordersize, pspc1, store, padding=[0.0, 0.5, 0.0])
      ! call closepacked(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, &
      !   bordersize, pspc1, store)
      ! call random([brdx1, brdx2, brdy1, brdy2, brdz1, brdz2], &
        ! bordersize, pspc1, store, displacement=.8)
      ! call createFixedBorders(store, ebc_x)
    case (ett_pulse, ett_ring)
      call setmhdmagneticpressure(1.)
      call setdiffisotropic("no")
      ! call setdiffisotropic("yes")
      call setdiffconductivity(1.)
      rho1 = 1.
      prs1 = 1.
      gamma = 5./3.

      brdx1 = -1.
      brdx2 =  1.
      if (resol == 0) then
        resol = int((brdx2-brdx1)/pspc1)
      end if
      pspc1 = (brdx2-brdx1)/resol
      pspc2 = pspc1
      brdy1 = -d2null*1.
      brdy2 =  d2null*1.
      brdz1 = -d3null*1.
      brdz2 =  d3null*1.
      bordersize = nb*pspc2
      ! call uniformV4(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, bordersize, pspc1, store, padding=0.5)
      ! call uniformV4(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, bordersize, pspc1, store, padding=0.5)
      call random([brdx1, brdx2, brdy1, brdy2, brdz1, brdz2], &
        bordersize, pspc1, store, displacement=.8)
      ! call closepacked(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, bordersize, pspc1, store)
      ! call createFixedBorders(store, ebc_all)
    case (ett_mti)
      call setmhdmagneticpressure(1.)
      call setdiffisotropic("yes")
      call setdiffconductivity(0.01)

      rho1  = 1.
      gamma = 5./3.

      brdx1 = 0.
      brdx2 = 1./10.
      if (resol == 0) then
        resol = int((brdx2-brdx1)/pspc1)
      end if
      pspc1 = (brdx2-brdx1)/resol
      pspc2 = pspc1
      brdy1 = 0.
      brdy2 = d2null*1./10.
      brdz1 = 0.
      brdz2 = d3null*1./10.
      bordersize = nb*pspc2
      call uniformV4(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, bordersize, pspc1, store, padding=0.5)
      call createFixedBorders(store, ebc_y)
      call getPartNumber(rpn,fpn)
      call set_density_profile(rpn+fpn,store,&
        brdy1-bordersize,brdy2+bordersize,&
        rhofunc=MTIHopkins2017,coord=2)
    case (ett_mtilowres)
      call setmhdmagneticpressure(1.)
      call setdiffisotropic("yes")
      ! call setdiffisotropic(0)
      call setdiffconductivity(0.1)

      rho1  = 1.
      gamma = 5./3.

      brdx1 = -1.
      ! brdx1 =  0.
      brdx2 =  1.
      if (resol == 0) then
        resol = int((brdx2-brdx1)/pspc1)
      end if
      pspc1 = (brdx2-brdx1)/resol
      pspc2 = pspc1
      ! brdy1 =  0.
      ! brdy2 =  1.
      brdy1 = -2.
      brdy2 =  2.
      brdz1 = 0.
      brdz2 = d3null*1.
      bordersize = nb*pspc2
      call uniformV4(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, bordersize, pspc1, store, padding=0.5)
      ! call createFixedBorders(store, ebc_y)
      call getPartNumber(rpn,fpn)
      ! call set_density_profile(rpn+fpn,store,&
      !   -2.,2.,&
      !   rhofunc=MTILowresHopkins2017,coord=2)
      ! call set_density_profile(rpn+fpn,store,&
      !   brdy1-bordersize,brdy2+bordersize,&
      !   rhofunc=MTIHopkins2017,coord=2)
    case (ett_soundwave)
      rho1 = 1.
      gamma= 5./3.
      eA = 0.1
      prs1 = 1.

      brdx1 = -1.
      brdx2 =  1.
      if (resol == 0) then
        resol = int((brdx2-brdx1)/pspc1)
      end if
      pspc1 = (brdx2-brdx1)/resol
      pspc2 = pspc1
      brdy1 = -pspc1*nb*2*d2null
      brdy2 =  pspc1*nb*2*d2null
      brdz1 = -pspc1*nb*2*d3null
      brdz2 =  pspc1*nb*2*d3null
      bordersize = nb*pspc2
      call uniformV4(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, bordersize, pspc1, store, padding=0.5)
    case (ett_hydroshock)
      call setmhdmagneticpressure(1.)
      call setdiffisotropic("yes")
      call setdiffconductivity(1.)

      gamma= 5./3.
      prs1 = 1.
      prs2 = prs1 / 10.
      rho1 = 1.
      rho2 = rho1 / 8.
      brdx1 = -.5
      brdx2 = .5
      if (resol == 0) then
        resol = abs(int(brdx1/pspc1))
      end if
      pspc1 = abs(brdx1/resol)
      pspc2 = pspc1*4.
      if (dim == 1) then
        pspc2 = pspc1 * 8.
      end if
      brdy1 = -d2null*pspc2*nb*2
      brdy2 =  d2null*pspc2*nb*2
      brdz1 = -d3null*pspc2*nb*2
      brdz2 =  d3null*pspc2*nb*2
      bordersize = nb*pspc2
      call uniformV4(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, bordersize, pspc1, store, dxmax=pspc2, padding=0.5)
      call createFixedBorders(store, ebc_x)
    case (ett_alfvenwave)
      call setmhdmagneticpressure(1.)
      call setdiffisotropic("no")
      call setdiffconductivity(1.)
      ! thetta
      theta = pi/2.
      phi   = pi/2.

      rho1  = 1.
      gamma = 5./3.
      eA    = 0.1
      prs1  = 0.1

      brdx1 = 0.
      ! brdx2 = 1./cos(theta)
      brdx2 = 1.
      if (resol == 0) then
        resol = int((brdx2-brdx1)/pspc1)
      end if
      pspc1 = (brdx2-brdx1)/resol
      pspc2 = pspc1
      nb = int(kr*hfac*2)
      if (dim > 1) then
        brdy1 = 0.
        brdy2 = 1./sin(theta)
      else
        brdy1 = 0.
        brdy2 = 0.
      end if
      if (dim == 3) then
        brdz1 = 0.
        brdz2 = 1.
      else
        brdz1 = 0.
        brdz2 = 0.
      end if
      bordersize = nb*pspc2
      call uniformV4(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, bordersize, pspc1, store, dxmax=pspc2, padding=0.5)
    case (ett_OTvortex)
      call setmhdmagneticpressure(1.)
      call setdiffisotropic("no")
      call setdiffconductivity(1.)

      prs1 = 0.133
      rho1 = 0.221
      gamma   = 5./3.

      brdx1 = 0.
      brdx2 = 1.
      if (resol == 0) then
        resol = int((brdx2-brdx1)/pspc1)
      end if
      pspc1 = (brdx2-brdx1)/resol
      pspc2 = pspc1
      nb = int(kr*hfac*2)
      brdy1 = 0. * d2null
      brdy2 = 1. * d2null
      brdz1 = -pspc1 * nb * d3null
      brdz2 =  pspc1 * nb * d3null
      bordersize = nb*pspc2
      call uniformV4(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, bordersize, pspc1, store, padding=0.5)
    case (ett_boilingtank)
      call setdiffisotropic("yes")
      call setdiffconductivity(1.)

      rho1  = 1.
      gamma = 5./3.

      brdx1 = -0.5
      brdx2 =  0.5
      if (resol == 0) then
        resol = int((brdx2-brdx1)/pspc1)
      end if
      pspc1 = (brdx2-brdx1)/resol
      pspc2 = pspc1
      brdy1 = 0.
      brdy2 = d2null*1.
      brdz1 = 0.
      brdz2 = d3null*1.
      bordersize = nb*pspc2
      call uniformV4(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, bordersize, pspc1, store, padding=0.5)
      call createFixedBorders(store, ebc_y)
      call getPartNumber(rpn,fpn)
      ! call set_density_profile(rpn+fpn,store,&
      !   brdy1-bordersize,brdy2+bordersize,&
      !   rhofunc=MTIHopkins2017,coord=2)
    case (ett_fld_gauss)
      include "setups/place_fld_gauss.f90"
    case default
      call error('Problem was not set', ivt, __FILE__, __LINE__)
    end select

    if (2*nb >= resol) then
      call error('There is not enough resolution to cover artificial particles', resol, __FILE__, __LINE__)
    end if

    call getPartNumber(rpn,fpn)
    call getPeriodPartNumber(ppn)
    call setGamma(gamma)

    lx = brdx2-brdx1
    ly = brdy2-brdy1
    lz = brdz2-brdz1

    if (ly == 0.) ly = 1.
    if (lz == 0.) lz = 1.

    n = rpn + fpn
    !$omp parallel do default(none)&
    !$omp private(i, sp, cca, qmatr)&
    !$omp private(ra)&
    !$omp shared(n, pspc1, pspc2, store, gamma, rho1, rho2)&
    !$omp shared(dim, hfac, tt, prs1, prs2, cv)&
    !$omp shared(brdx2, brdx1, brdy2, brdy1, brdz2, brdz1, lx, ly, lz)&
    !$omp shared(ivt, eA, kt, cs, theta, phi, d2null, d3null)&
    !$omp shared(rpn)
    do i=1,n
      ra(1) = store(es_rx,i)
      ra(2) = store(es_ry,i)
      ra(3) = store(es_rz,i)
      store(es_c,i) = cv
      store(es_om,i) = 1.
      sp = merge(pspc1, pspc2, ra(1) < 0)
      select case (ivt)
      case(ett_sin3)
        store(es_h,i) = hfac * sp
        store(es_m,i) = (sp**dim) * rho1
        store(es_den,i) = rho1
        ! isotropic case
        ! kcf(1,1,i) = kcf1
        ! kcf(2,2,i) = kcf1
        ! kcf(3,3,i) = kcf1
        ! if (pos(1,i) < 0) then
        !   kcf(i) = 1
        ! else
        !   kcf(i) = 10
        ! end if
        store(es_t,i) = 0
        if ( dim == 1) then
          store(es_t,i)  = sin(pi * (ra(1) - brdx1) / abs(brdx2-brdx1))
        elseif ( dim == 2 ) then
          store(es_t,i)  = sin(pi * (ra(1) - brdx1) / abs(brdx2-brdx1)) * &
                   sin(pi * (ra(2) - brdy1) / abs(brdy2-brdy1))
        elseif ( dim == 3 ) then
          store(es_t,i)  = sin(pi * (ra(1) - brdx1) / abs(brdx2-brdx1)) * &
                   sin(pi * (ra(2) - brdy1) / abs(brdy2-brdy1)) * &
                   sin(pi * (ra(3) - brdz1) / abs(brdz2-brdz1))
        end if
      case(ett_shock12)
        store(es_h,i) = hfac * sp
        ! store(es_m,i) = (sp**dim) * rho1
        store(es_m,i) = lx*ly*lz*rho1/rpn
        store(es_den,i) = rho1
        store(es_p,i) = prs1
        store(es_bx,i) = 1.
        store(es_by,i) = 0.
        store(es_bz,i) = 0.

        if (ra(1) < 0.) then
          store(es_t,i) = 1.
        else if (ra(1) > 0.) then
          store(es_t,i) = 2.
        else
          store(es_t,i) = 2.
          ! store(es_t,i) = 1.5
        end if
        ! store(es_u,i)  = store(es_t,i)
        store(es_u,i)  = prs1/(gamma -1)/rho1
      case(ett_pulse)
        store(es_h,i) = hfac * sp
        ! store(es_m,i) = (sp**dim) * rho1
        store(es_m,i) = (brdx2-brdx1)*(brdy2-brdy1)*rho1/rpn
        store(es_den,i) = rho1
        store(es_p,i) = prs1

        store(es_bx,i) = 1.
        store(es_by,i) = 0.
        store(es_bz,i) = 0.
        ! store(es_t,i)  = ((2.*pi)**(-dim/2.))/((0.1**2)**(dim/2.))*&
          ! exp(-0.5*(ra(1)*ra(1) + ra(2)*ra(2) + ra(3)*ra(3))/(0.1**2))
        ! store(es_u,i)  = store(es_t,i)
        store(es_kappa,i) = 1.
        !!!!!!!!!!!!!!!!!
        store(es_u,i)  = prs1/(gamma -1)/rho1
      case(ett_ring)
        store(es_h,i) = hfac * sp
        store(es_m,i) = (sp**dim) * rho1
        store(es_den,i) = rho1
        store(es_u,i)  = prs1/(gamma -1)/rho1
        store(es_p,i) = prs1

        ! initially in cylindric CS, but will be transphered in other system few lines lower
        store(es_bx,i) = 0.
        store(es_by,i) = 1.
        store(es_bz,i) = 0.

        cca(1) = sqrt(ra(1)*ra(1) + ra(2)*ra(2))
        cca(2) = atan(ra(2),ra(1))
        cca(3) = ra(3)

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

          store(es_bx:es_bz,i) = matmul(qmatr(:,:), store(es_bx:es_bz,i))
        end if
        ! exp[−(1/2)[(r−r0)^2/δr0^2 + φ^2/δφ0^2]]
        ! δr0 = 0.05 and r0 = 0.3 define a Gaussian ring
        !  at radius r0 of width δr, and δφ0 = 0.5
        ! is an initial Gaussian spread about φ=0intheφˆdirection
        store(es_t,i) = exp(-0.5*((cca(1) - 0.3)**2/(0.05*0.05) + cca(2)*cca(2)/(0.5*0.5)))
        store(es_u,i)  = store(es_t,i)
      case (ett_soundwave)
        store(es_h,i)   = hfac * sp
        store(es_den,i) = rho1 * (1. + eA * sin(pi * ra(1)))
        store(es_m,i)   = (sp**dim) * store(es_den,i)
        store(es_vx,i)  = eA*sin(pi*(ra(1)))
        store(es_p,i)   = prs1
        store(es_u,i)   = prs1/(gamma -1)/rho1
      case (ett_hydroshock)
        if (ra(1) <= 0.) then
          store(es_h,i)   = hfac * pspc1
          store(es_den,i) = rho1
          store(es_p,i)   = prs1
          store(es_m,i)   = (pspc1**dim) * rho1
          store(es_u,i)   = prs1/(gamma -1)/rho1
        else
          store(es_h,i)   = hfac * pspc2
          store(es_den,i) = rho2
          store(es_p,i)   = prs2
          store(es_m,i)   = (pspc2**dim) * rho2
          store(es_u,i)   = prs2/(gamma -1)/rho2
        end if
      case (ett_alfvenwave)
        store(es_h,i) = hfac * sp
        store(es_m,i) = (sp**dim) * rho1

        store(es_den,i) = rho1
        store(es_p,i)   = prs1
        store(es_u,i)   = prs1/(gamma -1)/rho1

        cca(:) = [cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta)]
        store(es_vx,i) = eA*0.
        store(es_vy,i) = eA*sin(2*pi*(dot_product(ra(:),cca(:))))
        store(es_vz,i) = eA*cos(2*pi*(dot_product(ra(:),cca(:))))
        store(es_bx:es_bz,i) = [1.,0.,0.] + store(es_vx:es_vz,i)
      case(ett_mti)
        store(es_h,i)   = hfac * sp
        store(es_m,i)   = (sp**dim) * rho1
        store(es_den,i) = rho1
        store(es_u,i)   = (3./2.)*(1. - ra(2)/3.)
        store(es_t,i)   = store(es_u,i)
        store(es_p,i)   = (gamma - 1.)*store(es_den,i)*store(es_u,i)
        store(es_bx:es_bz,i) = 0.
        store(es_vx:es_vz,i) = 0.
        cca(1) = dot_product(ra(:),ra(:))
        if (cca(1) > 0) then
          cca(:) = ra(:)/sqrt(cca(1))
          if (int(store(es_type,i)) == ept_real) then
            store(es_vy,i) = 1e-2*&
              sqrt(gamma * store(es_p,i) / store(es_den,i))*&
              sin(4.*pi*ra(1)/(brdx2-brdx1))
          end if
        end if
        store(es_bx,i) = 1e-11
      case(ett_mtilowres)
        store(es_h,i)   = hfac * sp
        store(es_m,i)   = (sp**dim) * rho1
        store(es_den,i) = rho1
        store(es_u,i)   = (3./2.)*(1. - abs(ra(2))/3.)
        store(es_t,i)   = store(es_u,i)
        store(es_p,i)   = (gamma - 1.)*store(es_den,i)*store(es_u,i)
        store(es_bx:es_bz,i) = 0.
        store(es_vx:es_vz,i) = 0.
        store(es_bx,i) = 1e-11
        ! if (int(store(es_type,i)) == ept_real) then
        !   cca(1) = dot_product(ra(:),ra(:))
        !   if (cca(1) > 0) then
        !     cca(:) = ra(:)/sqrt(cca(1))
        !       store(es_vy,i) = 1e-4*&
        !         sqrt(gamma * store(es_p,i) / store(es_den,i))*&
        !         sin(4.*pi*ra(1)/(brdx2-brdx1))
        !   end if
        ! end if
      case (ett_OTvortex)
        store(es_h,i)   = hfac * sp
        store(es_m,i)   = (sp**dim) * rho1
        store(es_den,i) = rho1
        store(es_p,i)   = prs1
        store(es_u,i)   = prs1/(gamma-1)/rho1

        store(es_vx,i) = -sin(2*pi*ra(2))
        store(es_vy,i) =  sin(2*pi*ra(1))
        store(es_vz,i) =  0.01
        store(es_bx,i) = -1./sqrt(4*pi)*sin(2*pi*ra(2))
        store(es_by,i) =  1./sqrt(4*pi)*sin(4*pi*ra(1))
        store(es_bz,i) =  0.
        if (dim == 1) then
          store(es_vy:es_vz,i) = 0.
          store(es_by:es_bz,i) = 0.
        elseif (dim == 2) then
          store(es_vz,i) = 0.
          store(es_bz,i) = 0.
        end if
      case (ett_boilingtank)
        store(es_h,i)   = hfac * sp
        store(es_m,i)   = (sp**dim) * rho1
        store(es_den,i) = rho1
        store(es_u,i)   = 1.
        store(es_t,i)   = 1.
        if (store(es_type,i) == ept_fixed) then
          if (ra(2) < brdy1) then
            if ((ra(1) > -0.075).and.(ra(1) < 0.075)) then
              store(es_t,i) = 10.
              store(es_u,i) = 10.
            end if
          end if
        end if
        store(es_p,i)   = (gamma - 1.)*store(es_den,i)*store(es_u,i)
        store(es_bx:es_by,i) = 1.
        store(es_vx:es_vz,i) = 0.
      case(ett_fld_gauss)
        include "setups/init_fld_gauss.f90"
      case default
        call error('There is no such initial condition', ivt, __FILE__, __LINE__)
      end select
      ! call printstore(store,i)
      ! read*
    end do
    !$omp end parallel do

    select case(ivt)
    case(ett_mti, ett_mtilowres)
      call setPartNumber(r=rpn+fpn,f=0)
      call clearPeriodicParticles(store)
      call findInsideBorderParticles(store)
      call createPeriodicBorder(store, ebc_x)
      call findneighboursKDT(store)
      call c1(store)
      call setPartNumber(r=rpn, f=fpn)
    end select

    call system_clock(finish)
    call addTime(' ic', finish - start)
  end subroutine setupV2

  ! subroutine printstore(s, i)
  !   real, allocatable, intent(in) :: s(:,:)
  !   integer, intent(in) :: i
  !   print*, "========"
  !   print*, "i  :", i, "type: ", s(es_type, i)
  !   print*, "pos: ", s(es_rx:es_rz, i)
  !   print*, "vel: ", s(es_vx:es_vz, i)
  !   print*, "acc: ", s(es_ax:es_az, i)
  !   print*, "P  : ", s(es_p, i)
  ! end subroutine printstore
end module
