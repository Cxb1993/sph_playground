module placer
  ! use omp_lib
  !
  ! use ArrayResize,      only: resize
  ! use kernel,           only: get_krad, &
  !                             get_w
  !                             getddwtype,&
  !                             getdim,&
  !                             ginitvar,&
  !                             gcoordsys,&
  !                             setdiffisotropic, &
  !                             setdiffconductivity, &
  !                             setmhdmagneticpressure,&
  !                             gethfac,&
  !                             getresolution,&
  !                             getspacing,&
  !                             getPartNumber,&
  !                             setPartNumber,&
  !                             setGamma
  ! use BC,               only: createFixedBorders, &
  !                             getPeriodPartNumber,&
  !                             createPeriodicBorder,&
  !                             findInsideBorderParticles,&
  !                             clearPeriodicParticles
  ! use neighboursearch,  only: getneighbours, &
  !                             getNeibListL1, &
  !                             getNeibListL2, &
  !                             findneighboursKDT
  ! use stretchmap,       only: set_density_profile
  ! use rhofuncs
  ! use circuit1, only: c1
  ! use circuit2, only: c2

  use const
  use state,            only: getStateVal,&
                              setStateVal
  use placeUniform,     only: uniformV4
  use placeClose,       only: closepacked
  use placeRandom,      only: random
  use kernel,           only: get_krad
  use errprinter,       only: error
  use timing,           only: addTime

  implicit none

  public :: place

  private

  integer(8) :: start=0, finish=0

contains

  subroutine place(store)
    real, allocatable, intent(inout) ::&
      store(:,:)
    integer ::&
      null2, null3, nb
    real ::&
      bordersize, kr, dim, placement, hfac,&
      xmin, xmax, ymin, ymax, zmin, zmax,&
      resolution, spacing, probbox(6)

    call system_clock(start)

    call getStateVal(ec_dim, dim)
    call getStateVal(ec_placement, placement)
    call getStateVal(ec_hfac, hfac)
    call getStateVal(ec_xmin, xmin)
    call getStateVal(ec_xmax, xmax)
    call getStateVal(ec_ymin, ymin)
    call getStateVal(ec_ymax, ymax)
    call getStateVal(ec_zmin, zmin)
    call getStateVal(ec_zmax, zmax)
    call getStateVal(ec_resolution, resolution)
    call getStateVal(ec_spacing, spacing)

    call get_krad(kr)

    null2 = 1
    null3 = 1
    if (int(dim) == 1) then
      null2 = 0
      null3 = 0
    else if (int(dim) == 2) then
      null3 = 0
    end if
    nb = int(kr*hfac)+1

    if ((int(resolution) <= 0).and.(spacing <= 0)) then
      call error('Both resolution and spacing were not specified', '', __FILE__, __LINE__)
    end if

    ! adjust spacing to follow resolution or to fit integer number of particles
    if (int(resolution) <= 0) then
      resolution = int((xmax-xmin)/spacing)
    end if
    spacing    = (xmax-xmin)/resolution
    bordersize = nb*spacing

    call setStateVal(ec_resolution, resolution)
    call setStateVal(ec_spacing, spacing)

    ymin = ymin * null2
    ymax = ymax * null2
    zmin = zmin * null3
    zmax = zmax * null3

    probbox(1:6) = [xmin, xmax, ymin, ymax, zmin, zmax]

    select case(int(placement))
    case(epl_uniform)
      call uniformV4(xmin, xmax, ymin, ymax, zmin, zmax, bordersize, spacing, store, padding=0.5)
    case(epl_random)
      call random(probbox, bordersize, spacing, store, displacement=.8)
    case(epl_closepacked)
      call closepacked(xmin, xmax, ymin, ymax, zmin, zmax, &
        bordersize, spacing, store, padding=[0.0, 0.5, 0.5])
    end select

    call system_clock(finish)
    call addTime('placing', finish - start)
  end subroutine place
end module
