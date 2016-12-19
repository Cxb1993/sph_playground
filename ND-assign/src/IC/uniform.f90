module uniform
  use kernel
  use utils
  use BC

  implicit none

  public :: make_uniform

  private

contains

  subroutine make_uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos)
    real, allocatable, intent(inout) :: pos(:,:)
    real, intent(in)     :: brdx1, brdx2, brdy1, brdy2, brdz1, brdz2
    real, intent(inout)  :: pspc1, pspc2
    integer, allocatable :: bX1(:), bY1(:), bZ1(:), bX2(:), bY2(:), bZ2(:)
    integer              :: dim, nb, bdx, bdy, bdz, ibx, iby, ibz, freeflag, freenumber, &
                            i, j, k, ix, iy, iz, n, nbnewX1, nbnewY1, nbnewZ1, nbnewX2, nbnewY2, nbnewZ2
    real                 :: spx, spy, spz, x, y, z

    call get_dim(dim)

    allocate(bX1(1))
    allocate(bX2(1))
    allocate(bY1(1))
    allocate(bY2(1))
    allocate(bZ1(1))
    allocate(bZ2(1))
    allocate(pos(3,1))

    n = 1
    nbnewX1 = 1
    nbnewY1 = 1
    nbnewZ1 = 1
    nbnewX2 = 1
    nbnewY2 = 1
    nbnewZ2 = 1


    ix = int((brdx2-brdx1)/pspc1)
    spx = merge(0.,(brdx2-brdx1)/ix, ix == 0)
    iy = int((brdy2-brdy1)/pspc1)
    spy = merge(0.,(brdy2-brdy1)/iy, iy == 0)
    iz = int((brdz2-brdz1)/pspc1)
    spz = merge(0.,(brdz2-brdz1)/iz, iz == 0)
    pspc1 = spx
    pspc2 = spx

    if (nb > 0) then
      bdx = nb
      bdy = merge(nb, 0, dim > 1)
      bdz = merge(nb, 0, dim == 3)
      ibx = 0
      iby = 0
      ibz = 0
    else
      bdx = 0
      bdy = 0
      bdz = 0
      ibx = abs(nb)
      iby = abs(nb)
      ibz = abs(nb)
    end if
    call set_sqare_box_sides(ix+1+2*bdx, iy+1+2*bdy, iz+1+2*bdz)
    freenumber = 0
    do i = (0-bdx),(ix+bdx)
      x = brdx1 + i * spx
      do j = (0-bdy),(iy+bdy)
        y = brdy1 + j * spy
        do k = (0-bdz),(iz+bdz)
          freeflag = 0
          z = brdz1 + k * spz
          if (x < brdx1 + ibx * spx) then
            freeflag = 1.
            if (size(bX1) < nbnewX1) then
              call resize(bX1,size(bX1),size(bX1)*2)
            end if
            bX1(nbnewX1) = n
            nbnewX1 = nbnewX1 + 1
          else if (x > brdx2 - ibx * spx) then
            freeflag = 1.
            if (size(bX2) < nbnewX2) then
              call resize(bX2,size(bX2),size(bX2)*2)
            end if
            bX2(nbnewX2) = n
            nbnewX2 = nbnewX2 + 1
          end if
          if (dim > 1) then
            if (y < brdy1 + iby * spy) then
              freeflag = 1.
              if (size(bY1) < nbnewY1) then
                call resize(bY1,size(bY1),size(bY1)*2)
              end if
              bY1(nbnewY1) = n
              nbnewY1 = nbnewY1 + 1
            else if (y > brdy2 - iby * spy) then
              freeflag = 1.
              if (size(bY2) < nbnewY2) then
                call resize(bY2,size(bY2),size(bY2)*2)
              end if
              bY2(nbnewY2) = n
              nbnewY2 = nbnewY2 + 1
            end if
            if (dim == 3) then
              if (z < brdz1 + ibz * spz) then
                freeflag = 1.
                if (size(bZ1) < nbnewZ1) then
                  call resize(bZ1,size(bZ1),size(bZ1)*2)
                end if
                bZ1(nbnewZ1) = n
                nbnewZ1 = nbnewZ1 + 1
              else if (z > brdz2 - ibz * spz) then
                freeflag = 1.
                if (size(bZ2) < nbnewZ2) then
                  call resize(bZ2,size(bZ2),size(bZ2)*2)
                end if
                bZ2(nbnewZ2) = n
                nbnewZ2 = nbnewZ2 + 1
              end if
            end if
          end if
          if (size(pos, dim=2) < n) then
            call resize3r(pos,size(pos, dim=2),size(pos, dim=2)*2)
          end if
          pos(1,n) = x
          pos(2,n) = y
          pos(3,n) = z
          n = n + 1
          freenumber = merge(freenumber + 1, freenumber, freeflag == 0)
        end do
      end do
    end do

    nbnewX1 = nbnewX1 - 1
    nbnewY1 = nbnewY1 - 1
    nbnewZ1 = nbnewZ1 - 1
    nbnewX2 = nbnewX2 - 1
    nbnewY2 = nbnewY2 - 1
    nbnewZ2 = nbnewZ2 - 1
    n = n - 1

    call resize3r(pos,n,n)

    write(*, "(A, F7.5, A, F7.5)") " # #   actual dx:   x1=", pspc1, "   x2=", pspc2
    print *, '# #      placed:', n
    print *, '# #  freenumber:', freenumber
    print *, '# #    border-x:', nbnewX1, nbnewX2
    print *, '# #    border-y:', nbnewY1, nbnewY2
    print *, '# #    border-z:', nbnewZ1, nbnewZ2

    call set_particles_numbers(n, abs(nb))
    call set_border(11, nbnewX1, bX1)
    call set_border(12, nbnewX2, bX2)
    call set_border(21, nbnewY1, bY1)
    call set_border(22, nbnewY2, bY2)
    call set_border(31, nbnewZ1, bZ1)
    call set_border(32, nbnewZ2, bZ2)

  end subroutine make_uniform
end module uniform
