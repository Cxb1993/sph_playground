module errcalc
  use errprinter,       only: error, warning
  use const
  use omp_lib
  use state,            only: getdim,&
                              get_equations,&
                              ginitvar,&
                              getdiffisotropic,&
                              getdiffconductivity,&
                              getPartNumber
  use neighboursearch,  only: getNeibListL1

  implicit none

  public :: err_sinxet,&
            err_diff_laplace, err_diff_graddiv, shockTube,&
            soundwaveperturbation_density, &
            soundwaveperturbation_velocity, &
            diff_artvisc, alfvenwave, hcpulse,&
            hcring, hcshock12

  private

contains

  subroutine shockTube(store, t, err)
    use exactshocktube
    real, allocatable, intent(in)    :: store(:,:)
    real, allocatable, intent(inout) :: err(:)
    real, intent(in)                 :: t

    integer           :: i, n
    real, allocatable :: exact(:), xpass(:)

    call getPartNumber(r=n)
    allocate(xpass(n))
    allocate(exact(n))
    xpass(:) = store(es_rx,1:n)
    call exact_shock(1, t, 1.4, 0., 1., 1./8., 1., 0.1, 0., 0., xpass, exact)

    do i = 1,n
      err(i) = abs(store(es_den,i)-exact(i))
    end do
  end subroutine

  subroutine hcpulse(store, t, err)
    real, allocatable, intent(in)    :: store(:,:)
    real, allocatable, intent(inout) :: err(:)
    real, intent(in)                 :: t

    integer, allocatable :: nlista(:)
    integer             :: i, j, dim, ivt, diso
    real :: &
      exact, x(3), num, dcnd, eps, pifac, ekt, ektd2, ektsq, epsdm1

    call getdim(dim)
    call getNeibListL1(nlista)
    call ginitvar(ivt)
    call getdiffisotropic(diso)
    call getdiffconductivity(dcnd)

    err(:) = 0.
    eps = 0.1
    pifac =(2*pi)**(-dim/2.)
    ekt = eps**2 + 2*dcnd*t
    ektd2 = ekt**(dim/2.)
    ektsq = ekt**(1./2.)
    epsdm1= eps**(dim-1)

    do j = 1,size(nlista)
      i = nlista(j)
      exact = 0.
      x(:) = store(es_rx:es_rz,i)
      num  = store(es_t,i)
      if (diso==1) then
        exact = pifac/ektd2*&
          exp(-0.5*((x(1)*x(1)+x(2)*x(2)+x(3)*x(3))/ekt))
      else
        exact = pifac/epsdm1/ektsq*&
          exp(-0.5*(x(1)*x(1)/ekt + (x(2)*x(2)+x(3)*x(3))/eps/eps))
          ! exp(-0.5*(x(2)*x(2)/ekt + (x(1)*x(1)+x(3)*x(3))/eps/eps))
      end if
      err(i) = abs(exact - num)
      ! print*, exact, num
      ! read*
      ! err(i) = dot_product(exact(1) - num(1,i), exact(1) - num(1,i))
    end do
  end subroutine hcpulse

  subroutine hcring(store, t, err)
    real, allocatable, intent(in)    :: store(:,:)
    real, allocatable, intent(inout) :: err(:)
    real, intent(in)                 :: t

    integer, allocatable :: nlista(:)
    integer             :: i, j, dim, ivt, diso
    real :: &
      exact, x(3), num, dcnd, eps, pifac,&
      ekt, ektd2, ektsq, epsdm1, rho, phi

    call getdim(dim)
    call getNeibListL1(nlista)
    call ginitvar(ivt)
    call getdiffisotropic(diso)
    call getdiffconductivity(dcnd)

    err(:) = 0.
    eps = 0.1
    pifac =(2*pi)**(-dim/2.)
    ekt = eps**2 + 2*dcnd*t
    ektd2 = ekt**(dim/2.)
    ektsq = ekt**(1./2.)
    epsdm1= eps**(dim-1)

    do j = 1,size(nlista)
      i = nlista(j)
      exact = 0.
      x(:) = store(es_rx:es_rz,i)
      rho = sqrt(x(1)*x(1) + x(2)*x(2))
      phi = atan(x(2),x(1))
      exact = exp(-0.5*((rho - 0.3)**2/(0.05*0.05) + (phi*phi)/(0.5*0.5 + 2.*1.*t/rho/rho)))
      num  = store(es_t,i)
      err(i) = abs(exact - num)
    end do
  end subroutine hcring

  subroutine hcshock12(store, t, err)
    real, allocatable, intent(in)    :: store(:,:)
    real, allocatable, intent(inout) :: err(:)
    real, intent(in)                 :: t
    integer, allocatable :: nlista(:)
    integer :: &
      i, j, dim, ivt, diso
    real :: &
      tr, tl, tt, x(3), exact, num

    call getNeibListL1(nlista)
    call getdiffisotropic(diso)

    tr = 2.
    tl = 1.

    do j = 1,size(nlista)
      i = nlista(j)
      exact = 0.
      x(:) = store(es_rx:es_rz,i)
      if (diso == 1) then
        tt = t
      else
        if (store(es_bx,i) == 1.) then
          tt = t
        else
          tt = 1e-15
        end if
      end if
      exact = (tr+tl)/2. + (tr-tl)/2. *erf(x(1)/sqrt(4*tt))
      num  = store(es_t,i)
      err(i) = abs(exact - num)
    end do
  end subroutine hcshock12

  subroutine err_sinxet(store, t, err)
    real, allocatable, intent(in)    :: store(:,:)
    real, allocatable, intent(inout) :: err(:)
    real, intent(in)                 :: t

    integer, allocatable :: nlista(:)
    integer             :: i, j, dim, ivt
    real                :: exact(3), x(3), num

    call getdim(dim)
    call getNeibListL1(nlista)
    call ginitvar(ivt)
    err(:) = 0.

    do j = 1,size(nlista)
      i = nlista(j)
      exact(:) = 0.
      x(:) = store(es_rx:es_rz,i)
      num  = store(es_t,i)
      if (ivt == 1) then
        if ( dim == 1 ) then
          exact(1) = sin(pi * (x(1) + 1.) / 2.) * exp(-(pi/2.)**2 * t)
        elseif ( dim == 2 ) then
          exact(1) = sin(pi * (x(1) + 1.) / 2.) * &
                  sin(pi * (x(2) + 1.) / 2.) * exp(-2 * (pi/2.)**2 * t)
        elseif ( dim == 3 ) then
          exact(1) = sin(pi * (x(1) + 1.) / 2.) * &
                  sin(pi * (x(2) + 1.) / 2.) * &
                  sin(pi * (x(3) + 1.) / 2.) * exp(-3 * (pi/2.)**2 * t)
        end if
        err(i) = abs(exact(1) - num)
      end if
      ! err(i) = dot_product(exact(1) - num(1,i), exact(1) - num(1,i))
    end do
  end subroutine err_sinxet

  subroutine alfvenwave(store, t, err)
    real, allocatable, intent(in)    :: store(:,:)
    real, allocatable, intent(inout) :: err(:)
    real, intent(in)                 :: t

    integer, allocatable :: nlista(:)
    integer              :: i, j, dim
    real                 :: exact, x1, B2, x(3), b(3)

    call getdim(dim)
    call getNeibListL1(nlista)
    err(:) = 0.
    do j = 1,size(nlista)
      i = nlista(j)
      x(:) = store(es_rx:es_rz,i)
      b(:) = store(es_bx:es_bz,i)
      x1 = (x(1) + 2*x(2) + 2*x(3))/3.
      B2 = (b(2) - 2*b(1))/sqrt(5.)

      exact = 0.1 * sin(2.*pi*(x1 - t))
      err(i) = abs(exact - B2)
    end do
  end subroutine alfvenwave

  subroutine soundwaveperturbation_density(store, t, err)
    real, allocatable, intent(in)    :: store(:,:)
    real, allocatable, intent(inout) :: err(:)
    real, intent(in)                 :: t

    integer, allocatable :: nlista(:)
    integer             :: i, j, dim
    real                :: exact

    call getdim(dim)
    call getNeibListL1(nlista)
    err(:) = 0.
    do j = 1,size(nlista)
      i = nlista(j)
      exact = 1. + 0.005 * sin(pi * (store(es_rx,i) - t))
      err(i) = abs(exact - store(es_den,i))
    end do
  end subroutine

  subroutine soundwaveperturbation_velocity(x, num, t, err)
    real, allocatable, intent(in)    :: num(:,:), x(:,:)
    real, allocatable, intent(inout) :: err(:)
    real, intent(in)                 :: t

    integer, allocatable :: nlista(:)
    integer             :: i, j, dim
    real                :: exact

    call getdim(dim)
    call getNeibListL1(nlista)
    err(:) = 0.
    !$omp parallel do default(none) &
    !$omp shared(x, num, err, dim, nlista, t) &
    !$omp private(exact, i, j)
    do j = 1,size(nlista)
      i = nlista(j)
      exact = 0.
      if ( dim == 1 ) then
        ! den
        ! exact = 1. + 0.005 * sin(pi * (x(1,i) - t))
        ! vel
        exact = 0.005 * sin(pi * (x(1,i) - t))
      end if
      err(i) = abs(exact - num(1,i))
      ! err(i) = dot_product(exact(1) - num(1,i), exact(1) - num(1,i))
    end do
    !$omp end parallel do
  end subroutine

  subroutine err_diff_laplace(x, num, err)
    real, allocatable, intent(in)    :: x(:,:), num(:,:)
    real, allocatable, intent(inout) :: err(:)

    integer, allocatable :: nlista(:)
    integer              :: i, j, dim
    real                 :: exact(1:3)

    call getdim(dim)
    call getNeibListL1(nlista)
    err(:) = 0.
    !$omp parallel do default(none) &
    !$omp shared(x, num, err, dim, nlista) &
    !$omp private(exact, i, j)
    do j = 1,size(nlista)
      i = nlista(j)
      exact(:) = 0.
      if ( dim == 1 ) then
        ! exact(1) = 2*Cos(x(1,i)) - x(1,i)*Sin(x(1,i))
        ! sin
        exact(1) = -sin(x(1,i))
      elseif ( dim == 2 ) then
        ! exact(1) = -x(2,i)*Sin(x(1,i))
        ! exact(2) = -x(1,i)*Sin(x(2,i))
        ! sin
        exact(1) = -sin(x(1,i))
        exact(2) = -sin(x(2,i))
      elseif ( dim == 3 ) then
        ! exact(1) = -(x(2,i)*Sin(x(1,i)))
        ! exact(2) = -(x(3,i)*Sin(x(2,i)))
        ! exact(3) = -(x(1,i)*Sin(x(3,i)))
        ! sin
        exact(1) = -sin(x(1,i))
        exact(2) = -sin(x(2,i))
        exact(3) = -sin(x(3,i))
      end if
      err(i) = dot_product(exact(:) - num(:,i),exact(:) - num(:,i))
    end do
    !$omp end parallel do
  end subroutine err_diff_laplace

  subroutine err_diff_graddiv(ptype, x, num, err)
    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(in)    :: x(:,:), num(:,:)
    real, allocatable, intent(inout) :: err(:)

    integer             :: n, i, dim, la
    real                :: exact(1:3)
    integer, allocatable :: nlista(:)

    call getdim(dim)
    n = size(ptype)
    err(:) = 0.

    call getNeibListL1(nlista)

    !$omp parallel do default(none) &
    !$omp shared(n,ptype, x,num,err,dim, nlista) &
    !$omp private(exact, i, la)
    do la = 1,size(nlista)
      i = nlista(la)
      ! print*, i, num(:,i)
      exact(:) = 0.
      if (dim == 1) then
        ! exact(1) = 0
        ! exact(1) = 2*Cos(x(1,i)) - (x(1,i))*Sin(x(1,i))
        ! sin
        exact(1) = -sin(x(1,i))
        ! grad only
        ! exact(1) = cos(x(1,i))
      end if
      if (dim == 2) then
        ! exact(1) = 1
        ! exact(2) = 1
        ! exact(1) = Cos(x(2,i)) - x(2,i)*Sin(x(1,i))
        ! exact(2) = Cos(x(1,i)) - x(1,i)*Sin(x(2,i))
        ! sin
        exact(1) = -sin(x(1,i))
        exact(2) = -sin(x(2,i))
        ! grad only
        ! exact(1) = cos(x(1,i))
        ! exact(2) = cos(x(2,i))
      end if
      if (dim == 3) then
        ! exact(1) = x(2,i) + x(3,i)
        ! exact(2) = x(1,i) + x(3,i)
        ! exact(3) = x(1,i) + x(2,i)
        ! exact(1) = Cos(x(3,i)) - (x(2,i)*Sin(x(1,i)))
        ! exact(2) = Cos(x(1,i)) - (x(3,i)*Sin(x(2,i)))
        ! exact(3) = Cos(x(2,i)) - (x(1,i)*Sin(x(3,i)))
        ! sin
        exact(1) = -sin(x(1,i))
        exact(2) = -sin(x(2,i))
        exact(3) = -sin(x(3,i))
        ! grad only
        ! exact(1) = cos(x(1,i))
        ! exact(2) = cos(x(2,i))
        ! exact(3) = cos(x(3,i))
      end if
      err(i) = dot_product(exact(:)-num(:,i),exact(:)-num(:,i))
    end do
    !$omp end parallel do
    call warning("bad norm specification. It is not L1", '', __FILE__, __LINE__)
  end subroutine err_diff_graddiv

  subroutine diff_artvisc(xin, num, err)
    real, allocatable, intent(in)    :: xin(:,:), num(:,:)
    real, allocatable, intent(inout) :: err(:)

    integer, allocatable :: nlista(:)
    integer              :: i, j, dim
    real                 :: exact(1:3), x(3)

    call getdim(dim)
    call getNeibListL1(nlista)
    err(:) = 0.
    !$omp parallel do default(none) &
    !$omp shared(xin, num, err, dim, nlista) &
    !$omp private(exact, i, j, x)
    do j = 1,size(nlista)
      i = nlista(j)
      x(:) = xin(:,i)
      exact(:) = 0.
      if ( dim == 1 ) then
        exact(1) = 2*Cos(x(1)) - x(1)*Sin(x(1)) + (2*Cos(x(1)) - x(1)*Sin(x(1)))/2.
        ! sin
        ! exact(1) = -3./2.*sin(x(1))
      elseif ( dim == 2 ) then
        exact(1) = 2*x(1)*Cos(x(2)) - (3*x(2)*Sin(x(1)))/2.
        exact(2) = Cos(x(1)) - x(1)**2*Sin(x(2)) + (2*Sin(x(2)) - x(1)**2*Sin(x(2)))/2.
        ! sin
        ! exact(1) = -3./2.*sin(x(1))
        ! exact(2) = -3./2.*sin(x(2))
      elseif ( dim == 3 ) then
        exact(1) = 3*x(1)**2*Cos(x(3)) - (3*x(2)*Sin(x(1)))/2.
        exact(2) = Cos(x(1)) - x(3)**2*Sin(x(2)) + (2*Sin(x(2)) - x(3)**2*Sin(x(2)))/2.
        exact(3) = 2*x(3)*Cos(x(2)) - x(1)**3*Sin(x(3)) + (6*x(1)*Sin(x(3)) - x(1)**3*Sin(x(3)))/2.
        ! sin
        ! exact(1) = -3./2.*sin(x(1))
        ! exact(2) = -3./2.*sin(x(2))
        ! exact(3) = -3./2.*sin(x(3))
      end if
      err(i) = dot_product(exact(:) - num(:,i),exact(:) - num(:,i))
      ! print*, err(i)
      ! print*, num(:,i)
      ! print*, exact
      ! read*
    end do
    !$omp end parallel do
    call warning("bad norm specification. It is not L1", '', __FILE__, __LINE__)
  end subroutine
end module
