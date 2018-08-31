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
    real, optional, intent(in)  :: padding

    integer              :: dim, n, ptsz, &
                             ix1, ix2, iy1, iy2, iz1, iz2, i, j, k, &
                             d2null, d3null
    real :: &
      spx, spy, spz, dmx, pdg

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
      pdg = 0.
    else
      pdg = padding
    end if

    ix1 = int(xmin/dx)
    ix2 = int(xmax/dx)
    do i = ix1, ix2 - int(2*pdg)
      spx = dx
      spy = sqrt(3./4.*spx*spx)
      spz = dx/2.
      iy1 = int(ymin/spy)
      iy2 = int(ymax/spy)
      iz1 = int(zmin/spy)
      iz2 = int(zmax/spy)
      do j = iy1, iy2 - int(2*pdg)*d2null
        do k = iz1, iz2 - int(2*pdg)*d3null
          if (.not.((i == ix1).and.(mod(j,2)==0))) then
            ptsz = size(store, dim=2)
            if (ptsz < n) then
              call resize(store, ptsz, ptsz*2)
            end if
            store(:,n) = 0.
            store(es_rx,n)    = spx*i + mod(j,2)*spz
            store(es_ry,n)    = spy*j
            store(es_rz,n)    = spy*k
            store(es_type,n)  = ept_real
            n = n + 1
          end if
        end do
      end do
    end do
    n = n - 1
    call setPartNumber(r=n,f=0)
    ! call setBorders(xmin, xmax, ymin, ymax, zmin, zmax, dbsz, pdg*sp)
  end subroutine closepacked
end module placeClose
