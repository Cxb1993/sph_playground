module placeRandom
  use const
  use kernel
  use state,        only: getdim, setBorders, setPartNumber
  use ArrayResize,  only: resize

  implicit none

  public :: random

  private

contains
  subroutine random( &
      boundBox, dbsz, dxmin, store, dxmax, padding, displacement)
    real, allocatable, intent(inout)    :: store(:,:)
    real, intent(in) :: &
      boundBox(6), dbsz
    real, intent(inout)         :: dxmin
    real, optional, intent(in)  :: dxmax, padding, displacement
    integer :: &
      dim, n, ptsz, &
      ix1, ix2, iy1, iy2, iz1, iz2, i, j, k, &
      d2null, d3null, resol
    real :: &
      sp, dmx, pdg, rnddr(3), drdis, &
      xmin, xmax, ymin, ymax, zmin, zmax

    call getdim(dim)
    d2null = 1
    d3null = 1
    if (dim == 1) then
      d2null = 0
      d3null = 0
    else if (dim == 2) then
      d3null = 0
    end if

    allocate(store(es_total,1))
    n = 1

    if (.not.present(dxmax)) then
      dmx = dxmin
    else
      dmx = dxmax
    end if
    if (.not.present(padding)) then
      pdg = 0.
    else
      pdg = padding
    end if
    if (.not.present(displacement)) then
      drdis = 0.
    else
      drdis = displacement
    end if

    xmin = boundBox(1)
    xmax = boundBox(2)
    ymin = boundBox(3)
    ymax = boundBox(4)
    zmin = boundBox(5)
    zmax = boundBox(6)

    ! call random_seed(put=18081991)

    ix1 = int(xmin/dxmin)
    ix2 = int(xmax/dmx)
    resol = ix2 - ix1
    do i = ix1, ix2 - int(2*pdg)
      if (i <= 0) then
        sp = dxmin
      else
        sp = dmx
      end if
      iy1 = int(ymin/sp)
      iy2 = int(ymax/sp)
      iz1 = int(zmin/sp)
      iz2 = int(zmax/sp)
      do j = iy1, iy2 - int(2*pdg)*d2null
        do k = iz1, iz2 - int(2*pdg)*d3null
          ptsz = size(store, dim=2)
          if (ptsz < n) then
            call resize(store, ptsz, ptsz*2)
          end if
          store(:,n) = 0.
          call random_number(rnddr)
          store(es_rx,n)    =        (pdg + i + drdis*rnddr(1))*sp
          store(es_ry,n)    = d2null*(pdg + j + drdis*rnddr(2))*sp
          store(es_rz,n)    = d3null*(pdg + k + drdis*rnddr(3))*sp

          if (store(es_rx,n)<xmin) then
            store(es_rx,n) = xmin + (xmax - xmin)/resol/4.
          end if
          if (store(es_rx,n)>xmax) then
            store(es_rx,n) = xmax - (xmax - xmin)/resol/4.
          end if
          if (store(es_ry,n)<ymin) then
            store(es_ry,n) = ymin + (ymax - ymin)/resol/4.
          end if
          if (store(es_ry,n)>ymax) then
            store(es_ry,n) = ymax - (ymax - ymin)/resol/4.
          end if
          if (store(es_rz,n)<zmin) then
            store(es_rz,n) = zmin + (zmax - zmin)/resol/4.
          end if
          if (store(es_rz,n)>zmax) then
            store(es_rz,n) = zmax - (zmax - zmin)/resol/4.
          end if
          store(es_type,n)  = ept_real
          n = n + 1
        end do
      end do
    end do
    n = n - 1
    call setPartNumber(r=n,f=0)
    call setBorders(xmin, xmax, ymin, ymax, zmin, zmax, dbsz, pdg*sp)
  end subroutine random
end module placeRandom
