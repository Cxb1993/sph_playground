module IC
  use omp_lib

  use timing,           only: addTime
  use ArrayResize,      only: resize
  use const
  use kernel,           only: get_krad, &
                              get_w
  use state,            only: get_tasktype,&
                              getkerntype,&
                              getdim,&
                              ginitvar,&
                              gcoordsys,&
                              sorigin, &
                              diff_conductivity,&
                              diff_isotropic, &
                              mhd_magneticconstant
  use BC,               only: realpartnumb, &
                              artpartnumb, &
                              bordersize, &
                              xmin, xmax, &
                              ymin, ymax, &
                              zmin, zmax
  use initpositions,    only: uniform, uniformV3
  use neighboursearch,  only: getneighbours, &
                              getNeibListL1, &
                              getNeibListL2, &
                              findneighboursKDT
  implicit none

  public :: setupIC

  private
  integer(8) :: start=0, finish=0
contains

  subroutine setupIC(n, sk, g, pspc1, resol, &
    pos, v, dv, mas, den, sln, prs, iu, du, cf, kcf, dcf, ptype)
    integer, allocatable, intent(inout) :: ptype(:)
    real, allocatable, intent(inout), dimension(:,:)  :: pos, v, dv, cf, dcf
    real, allocatable, intent(inout), dimension(:,:,:):: kcf
    real, allocatable, intent(inout), dimension(:)    :: mas, den, sln, prs, iu, du
    real, intent(in)        :: sk
    real, intent(inout)     :: pspc1, g
    integer, intent(inout)  :: n, resol

    real                 :: kr, prs1, prs2, rho1, rho2, kcf1, sp, &
                            brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, period, eA, &
                            cca(3), qmatr(3,3), qtmatr(3,3), pspc2, theta, phi
    integer              :: i, nb, tt, kt, dim, ivt, cs, d2null, d3null
    ! integer, allocatable :: nlista(:)

    call system_clock(start)
    ! do i=1,5
    !   allocate(pos(3,10**i))
    !   call system_clock(start)
    !   call resize(pos,10**i,10**(i*2))
    !   call system_clock(finish, count_rate=cr)
    !   deallocate(pos)
    !   print*, "empty  resize time: ", 10**i, "to ", 10**(i*2), (finish - start)/real(cr)
    !
    !   allocate(pos(3,10**i))
    !   CALL random_number(pos)
    !   call system_clock(start)
    !   call resize(pos,10**i,10**(i*2))
    !   call system_clock(finish, count_rate=cr)
    !   deallocate(pos)
    !   print*, "random resize time: ", 10**i, "to ", 10**(i*2), (finish - start)/real(cr)
    !
    !   allocate(pos(3,10**(i*2)))
    !   CALL random_number(pos)
    !   call system_clock(start)
    !   call resize(pos,10**(i*2), 10**i)
    !   call system_clock(finish, count_rate=cr)
    !   deallocate(pos)
    !   print*, "back   resize time: ", 10**(i*2), "to ", 10**i, (finish - start)/real(cr)
    !   print*, '====='
    ! end do
    ! read*
    call getkerntype(kt)
    call get_tasktype(tt)
    call get_krad(kr)
    call getdim(dim)
    call ginitvar(ivt)
    call gcoordsys(cs)

    d2null = 1
    d3null = 1
    if (dim == 1) then
      d2null = 0
      d3null = 0
    else if (dim == 2) then
      d3null = 0
    end if

    ! don't need to have it different for 2nw, just copy it periodicly
    nb = int(kr*sk) + 1

    select case(ivt)
    case(ett_sin3, ett_shock12, ett_pulse)
      print*, "Not set yet. FIX ME. IC.f90. line 103."
      stop
    case (ett_mti)
      mhd_magneticconstant = 0.00000125663706
      diff_isotropic = 1.
      diff_conductivity = 0.01

      rho1 = 1.
      g    = 5./3.

      brdx1 = 0.
      brdx2 = 1./10.
      if (pspc1 /= 0) then
        resol = int((brdx2-brdx1)/pspc1)
      end if
      pspc1 = (brdx2-brdx1)/resol
      pspc2 = pspc1
      nb = int(kr*sk*2)
      if (dim > 1) then
        brdy1 = 0.
        brdy2 = 1./10.
      else
        brdy1 = 0.
        brdy2 = 0.
      end if
      if (dim == 3) then
        brdz1 = 0.
        brdz2 = 1./10.
      else
        brdz1 = 0.
        brdz2 = 0.
      end if
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos, ptype, randomise=0.)
    case (ett_ring)
      diff_conductivity = 1.

      rho1 = 1.
      kcf1 = 1.

      brdx1 = -1.
      brdx2 =  1.
      if (pspc1 /= 0) then
        resol = int((brdx2-brdx1)/pspc1)
      end if
      pspc1 = (brdx2-brdx1)/resol
      pspc2 = pspc1
      if (dim > 1) then
          brdy1 = -1.
          brdy2 =  1.
      else
        brdy1 = 0.
        brdy2 = 0.
      end if
      if (dim == 3) then
          brdz1 = -1.
          brdz2 =  1.
      else
        brdz1 = 0.
        brdz2 = 0.
      end if
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos, ptype, randomise=0.)
    case (ett_soundwave)
      rho1 = 1.
      g = 5./3.
      eA = 0.1
      prs1 = 1.

      brdx1 = -1.
      brdx2 =  1.
      if (pspc1 /= 0) then
        resol = int((brdx2-brdx1)/pspc1)
      end if
      pspc1 = (brdx2-brdx1)/resol
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
    case (ett_hydroshock)
      g = 1.4
      prs1 = 1.
      prs2 = prs1 / 10.
      rho1 = 1.
      rho2 = rho1 / 8.
      if (pspc1 /= 0) then
        resol = int(brdx2/pspc1)
      end if
      pspc1 = brdx2/resol
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
    case (ett_alfvenwave)
      mhd_magneticconstant = 1.
      ! thetta
      theta = pi/2.
      phi = pi/2.

      rho1 = 1.
      g    = 5./3.
      eA   = 0.1
      prs1 = 0.1

      brdx1 = 0.
      ! brdx2 = 1./cos(theta)
      brdx2 = 1.
      if (pspc1 /= 0) then
        resol = int((brdx2-brdx1)/pspc1)
      end if
      pspc1 = (brdx2-brdx1)/resol
      pspc2 = pspc1
      nb = int(kr*sk*2)
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
      call uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos, ptype, randomise=0.)
    case (ett_OTvortex)
      mhd_magneticconstant = 1.
      diff_isotropic = 1.
      diff_conductivity = 0.

      prs1 =  5./(12.*pi)
      rho1 = 25./(36.*pi)
      g    = 5./3.

      brdx1 = -1
      brdx2 =  1
      if (resol == 0) then
        resol = int((brdx2-brdx1)/pspc1)
      end if
      pspc1 = (brdx2-brdx1)/resol
      pspc2 = pspc1
      nb = int(kr*sk*2)
      brdy1 = -1 * d2null
      brdy2 =  1 * d2null
      brdz1 = -pspc1 * nb * d3null
      brdz2 =  pspc1 * nb * d3null
      call uniformV3(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pos, padding=0.5)
      artpartnumb = (2*nb)*((2*resol+2*nb)*d2null+1)*((resol+2*nb)*d3null+1)
    case default
      print *, 'Problem was not set in IC.f90: line 220.'
      stop
    end select

    if (2*nb >= resol) then
      print*, '   +-> not enough resolution to cover artificial particles'
      stop
    end if

    n = realpartnumb+artpartnumb
    call resize(pos, realpartnumb, n)
    if (n < 2) then
      error stop 'There is only ' // char(n+48) // ' particles in the domain'
    else
      bordersize = nb*pspc1
      xmin = brdx1
      xmax = brdx2
      ymin = brdy1
      ymax = brdy2
      zmin = brdz1
      zmax = brdz2
      if (.not.(allocated(ptype))) then
        allocate(ptype(n))
        ptype(1:realpartnumb) = 1
        pos(:,realpartnumb+1:n) = 0.
        ptype(realpartnumb+1:n) = 0
      end if
      artpartnumb = 0
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

    !-------------------------------!
    !        particles values       !
    !-------------------------------!

    ! sln(:) = sk * pspc1
    ! call findneighboursKDT(ptype, pos, sln)

    !$omp parallel do default(none)&
    !$omp private(i, sp, cca, qmatr, qtmatr)&
    !$omp shared(n, pspc1, pspc2, pos, sln, den, prs, mas, iu, g, rho1, rho2)&
    !$omp shared(cf, kcf, dim, sk, tt, prs1, prs2, kcf1)&
    !$omp shared(brdx2, brdx1, brdy2, brdy1, brdz2, brdz1, v, period, ptype)&
    !$omp shared(ivt, eA, kt, cs, theta, phi)
    do i=1,realpartnumb
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
      ! if (ptype(i) /= 0) then
        sp = merge(pspc1, pspc2, pos(1,i) < 0)

        select case (ivt)
        case(ett_sin3)
          sln(i) = sk * sp
          mas(i) = (sp**dim) * rho1
          den(i) = rho1

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
        case(ett_shock12)
          sln(i) = sk * sp
          mas(i) = (sp**dim) * rho1
          den(i) = rho1

          kcf(1,1,i) = 0.
          kcf(1,2,i) = 1.
          kcf(1,3,i) = 0.
          ! cf(:, i) = sin(pos(:,i))
          if (pos(1,i) < 0.) then
            cf(1,i) = 1.
          else if (pos(1,i) > 0.) then
            cf(1,i) = 2.
          else
            cf(1,i) = 1.5
          end if
        case(ett_pulse)
          sln(i) = sk * sp
          mas(i) = (sp**dim) * rho1
          den(i) = rho1

          kcf(1,1,i) = 1.
          kcf(1,2,i) = 0.
          kcf(1,3,i) = 0.
          ! cf(:, i) = sin(pos(:,i))
          cf(1,i) = (2.*pi)**(-dim/2.)/(0.1**2)**(dim/2.)*&
                    exp(-0.5*(pos(1,i)*pos(1,i) + pos(2,i)*pos(2,i) + pos(3,i)*pos(3,i))/(0.1**2))
        case(ett_ring)
          sln(i) = sk * sp
          mas(i) = (sp**dim) * rho1
          den(i) = rho1
          ! initially in cylindric CS, but will be transphered in other system few lines lower
          kcf(1,1,i) = 0.
          kcf(2,1,i) = 1.
          kcf(3,1,i) = 0.

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

            kcf(:,1,i) = matmul(qmatr(:,:), kcf(:,1,i))
          end if
          ! exp[−(1/2)[(r−r0)^2/δr0^2 + φ^2/δφ0^2]]
          ! δr0 = 0.05 and r0 = 0.3 define a Gaussian ring
          !  at radius r0 of width δr, and δφ0 = 0.5
          ! is an initial Gaussian spread about φ=0intheφˆdirection
          cf(1,i) = exp(-0.5*((cca(1) - 0.3)**2/(0.05*0.05) + cca(2)*cca(2)/(0.5*0.5)))
        case (ett_soundwave)
          sln(i) = sk * sp
          den(i) = rho1 * (1. + eA * sin(pi * pos(1,i)))
          mas(i) = (sp**dim) * den(i)
          v(1,i) = eA*sin(pi*(pos(1,i)))
          prs(i) = prs1
          iu(i)  = prs1/(g-1)/rho1
        case (ett_hydroshock)
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
        case (ett_alfvenwave)
          sln(i) = sk * sp
          mas(i) = (sp**dim) * rho1

          den(i) = rho1
          prs(i) = prs1
          iu(i)  = prs1/(g-1)/rho1

          cca(:) = [cos(theta)*cos(phi), cos(theta)*sin(phi), sin(theta)]

          ! kcf(1,1,i) = 1.*cos(cca(1)) - cca(3) * sin(cca(1))
          ! kcf(2,1,i) = 1.*sin(cca(1)) + cca(3) * cos(cca(1))
          ! kcf(3,1,i) = eA*cos(2*pi*cca(2))
          ! v(1,i) = kcf(1,1,i)
          ! v(2,i) = kcf(2,1,i)
          ! v(3,i) = kcf(3,1,i)
          v(:,i) = eA*[0.,&
            sin(2*pi*(dot_product(pos(:,i),cca(:)))),&
            cos(2*pi*(dot_product(pos(:,i),cca(:))))]
          kcf(:,1,i) = [1.,0.,0.] + v(:,i)
          ! v(1,i) = 0.
          ! v(2,i) = kcf(2,1,i)
          ! v(3,i) = kcf(3,1,i)
          ! if (dim == 1) then
          !   v(2:3,i)   = 0.
          !   kcf(2:3,1,i)  = 0.
          ! elseif (dim == 2) then
          !   v(3,i)    = 0.
          !   kcf(3,1,i)  = 0.
          ! end if
        case(ett_mti)
          sln(i) = sk * sp
          mas(i) = (sp**dim) * rho1
          den(i) = rho1
          prs(i) = prs1
          iu(i)  = (3./2.)*(1. - pos(2,i)/3.)
          cca(1) = dot_product(pos(:,i),pos(:,i))
          if (cca(1) /= 0) then
            cca(:) = pos(:,i)/sqrt(cca(1))
            kcf(:,1,i) = 10e-11*cca(1)
            v(:,i) = 10e-2*1.*sin(4.*pi*pos(1,i)/(brdx2-brdx1))*cca(2)
          end if
          if (dim == 1) then
            v(2:3,i)   = 0.
            kcf(2:3,1,i)  = 0.
          elseif (dim == 2) then
            v(3,i)    = 0.
            kcf(3,1,i)  = 0.
          end if
        case (ett_OTvortex)
          sln(i) = sk * sp
          mas(i) = (sp**dim) * rho1
          den(i) = rho1
          prs(i) = prs1
          iu(i)  = prs1/(g-1)/rho1

          v(1,i) = -sin(2*pi*(pos(2,i)-brdy1))
          v(2,i) =  sin(2*pi*(pos(1,i)-brdx1))
          v(3,i) = 0.01
          kcf(1,1,i) = -1./sqrt(4*pi)*sin(2*pi*(pos(2,i)-brdy1))
          kcf(2,1,i) =  1./sqrt(4*pi)*sin(4*pi*(pos(1,i)-brdx1))
          kcf(3,1,i) =  0.
          if (dim == 1) then
            v(2:3,i)   = 0.
            kcf(2:3,1,i)  = 0.
          elseif (dim == 2) then
            v(3,i)    = 0.
            kcf(3,1,i)  = 0.
          end if
        end select
      ! end if
    end do
    !$omp end parallel do
    call system_clock(finish)
    call addTime(' ic', finish - start)
  end subroutine
end module
