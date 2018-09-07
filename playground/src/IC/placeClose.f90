module placeClose
  use const
  use kernel
  use state,        only: getdim, setBorders, setPartNumber
  use ArrayResize,  only: resize

  implicit none

  public :: closepacked

  private

contains
  subroutine closepacked(xmin, xmax, ymin, ymax, zmin, zmax, dbsz, dx, store, padding)
    real, allocatable, intent(inout)    :: store(:,:)
    real, intent(in) :: &
      xmin, xmax, ymin, ymax, zmin, zmax, dbsz
    real, intent(inout)         :: dx
    real, optional, intent(in)  :: padding(3)

    integer :: &
      dim, n, ptsz, ix1, ix2, iy1, iy2, iz1, iz2, i, j, k, &
      d2null, d3null
    real :: &
      spdx, sph, sph3, spdx2, dmx, pdg(3), &
      dymineps, dymaxeps, dzmineps, dzmaxeps, &
      newymin, newymax, newzmin, newzmax


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

    if (.not.present(padding)) then
      pdg(:) = 0.
    else
      pdg(:) = padding(:)
    end if


    ! spdx = merge(...)
    spdx = dx
    ix1 = int(xmin/spdx)
    ix2 = int(xmax/spdx)
    sph = sqrt(3./4.*spdx*spdx)
    spdx2 = spdx/2.
    sph3 = sph/3.

    iy1 = int(ymin/sph)
    iy2 = int(ymax/sph)
    dymineps = 0.
    dymaxeps = 1.
    newymin = sph*iy1
    newymax = sph*iy2
    if ((iy1 - iy2) /= 0) then
      dymineps = abs((ymax - ymin) - (iy2 - iy1)*sph)/abs(ymax - ymin)
      dymaxeps = abs((ymax - ymin) - (iy2 - iy1 + 2)*sph)/abs(ymax - ymin)
      if (dymaxeps < dymineps) then
        iy1  = iy1 - 1
        newymin = iy1 * sph
        iy2  = iy2 + 1
        newymax = iy2 * sph
      end if
    end if

    iz1 = int(zmin/sph)
    iz2 = int(zmax/sph)
    dzmineps = 0.
    dzmaxeps = 1.
    newzmin = zmin
    newzmax = zmax
    if ((iz1-iz2)/= 0) then
      dzmineps = abs((zmax - zmin) - (iz2 - iz1)*sph)/abs(zmax - zmin)
      dzmaxeps = abs((zmax - zmin) - (iz2 - iz1 + 2)*sph)/abs(zmax - zmin)
      if (dzmaxeps < dzmineps) then
        iz1  = iz1 - 1
        newzmin = iz1 * sph
        iz2  = iz2 + 1
        newzmax = iz2 * sph
      end if
    end if
    do i = ix1, ix2 - int(2*pdg(1))
      do j = iy1, iy2 - int(2*pdg(2))*d2null
        do k = iz1, iz2 - int(2*pdg(3))*d3null
          if (mod(k-iz1,2) == 0) then
            if (.not.((i == ix2).and.(mod(j-iy1,2)/=0))) then
              ptsz = size(store, dim=2)
              if (ptsz < n) then
                call resize(store, ptsz, ptsz*2)
              end if
              store(:,n) = 0.
              store(es_rx,n)    =        pdg(1)*spdx + spdx*i + mod(j-iy1,2)*spdx2
              store(es_ry,n)    = d2null*pdg(2)*sph +  sph*j
              store(es_rz,n)    = d3null*pdg(3)*sph +  sph*k
              store(es_type,n)  = ept_real
              n = n + 1
            end if
          else
            if (.not.((i == ix1).and.(mod(j-iy1,2)==0))) then
              ptsz = size(store, dim=2)
              if (ptsz < n) then
                call resize(store, ptsz, ptsz*2)
              end if
              store(:,n) = 0.
              store(es_rx,n)    = spdx*i - spdx2 + mod(j-iy1,2)*spdx2
              store(es_ry,n)    =  sph*j + sph3
              store(es_rz,n)    =  sph*k
              store(es_type,n)  = ept_real
              n = n + 1
            end if
          end if
        end do
      end do
    end do
    n = n - 1
    call setPartNumber(r=n,f=0)
    call setBorders(xmin, xmax, newymin, newymax, newzmin, newzmax, dbsz, pdg(2)*sph)
  end subroutine closepacked
end module placeClose
