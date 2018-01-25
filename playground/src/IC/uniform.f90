module uniform
  use const
  use kernel
  use state,        only: getdim, setBorders, setPartNumber
  use ArrayResize,  only: resize

  implicit none

  public :: uniformV4

  private

contains
  subroutine uniformV4(xmin, xmax, ymin, ymax, zmin, zmax, dbsz, dxmin, store, dxmax, padding)
    real, allocatable, intent(inout)    :: store(:,:)
    real, intent(in) :: &
      xmin, xmax, ymin, ymax, zmin, zmax, dbsz
    real, intent(inout)         :: dxmin
    real, optional, intent(in)  :: dxmax, padding

    integer              :: dim, n, ptsz, &
                             ix1, ix2, iy1, iy2, iz1, iz2, i, j, k, &
                             d2null, d3null
    real                 :: sp, dmx, pdg

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

    ix1 = int(xmin/dxmin)
    ix2 = int(xmax/dmx)
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
          store(es_rx,n)    =        pdg*sp + sp*i
          store(es_ry,n)    = d2null*pdg*sp + sp*j
          store(es_rz,n)    = d3null*pdg*sp + sp*k
          store(es_type,n)  = ept_real
          n = n + 1
        end do
      end do
    end do
    n = n - 1
    call setPartNumber(r=n,f=0)
    call setBorders(xmin, xmax, ymin, ymax, zmin, zmax, dbsz)
  end subroutine uniformV4
end module uniform
