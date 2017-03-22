module initpositions
  use kernel
  use utils
  use BC

  implicit none

  public :: place_uniform, place_close_packed_fcc

  private

contains

  subroutine place_uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos, ptype)
    real,    allocatable, intent(inout) :: pos(:,:)
    integer, allocatable, intent(inout) :: ptype(:)
    real, intent(inout)  :: pspc1, pspc2, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2

    integer, allocatable :: bX1(:), bY1(:), bZ1(:), bX2(:), bY2(:), bZ2(:)
    integer              :: dim, nb, bdx, bdy, bdz, freeflag, freenumber, &
                            i, j, k, ix, iy, iz, n, nbnewX1, nbnewY1, nbnewZ1, nbnewX2, nbnewY2, nbnewZ2,&
                            ptsz !particle array size
                            ! , ibx, iby, ibz ! number thet were used for 'niternal' border setting
    real                 :: spx, spy, spz, x, y, z, eps

    eps = 10*epsilon(0.)

    call get_dim(dim)

    allocate(bX1(1))
    allocate(bX2(1))
    allocate(bY1(1))
    allocate(bY2(1))
    allocate(bZ1(1))
    allocate(bZ2(1))
    allocate(pos(3,1))
    allocate(ptype(1))

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

    ! if (nb > 0) then
    bdx = nb
    bdy = merge(nb, 0, dim > 1)
    bdz = merge(nb, 0, dim == 3)
    !   ibx = 0
    !   iby = 0
    !   ibz = 0
    ! else
    !   bdx = 0
    !   bdy = 0
    !   bdz = 0
    !   ibx = abs(nb)
    !   iby = abs(nb)
    !   ibz = abs(nb)
    ! end if
    ! print *, bdx, ix, brdx1, brdx2, ibx
    call set_sqare_box_sides(ix+1+2*bdx, iy+1+2*bdy, iz+1+2*bdz)
    freenumber = 0
    do i = (0-bdx),(ix+bdx)
      x = brdx1 + i * spx
      do j = (0-bdy),(iy+bdy)
        y = brdy1 + j * spy
        do k = (0-bdz),(iz+bdz)
          freeflag = 0
          z = brdz1 + k * spz
          ! if (x < brdx1 + ibx * spx - eps) then
          if (x < brdx1 - eps) then
            freeflag = 1.
            if (size(bX1) < nbnewX1) then
              call resize(bX1,size(bX1),size(bX1)*2)
            end if
            bX1(nbnewX1) = n
            nbnewX1 = nbnewX1 + 1
          ! else if (x > brdx2 - ibx * spx + eps) then
          else if (x > brdx2 + eps) then
            freeflag = 1.
            if (size(bX2) < nbnewX2) then
              call resize(bX2,size(bX2),size(bX2)*2)
            end if
            bX2(nbnewX2) = n
            nbnewX2 = nbnewX2 + 1
          end if
          if (dim > 1) then
            ! if (y < brdy1 + iby * spy - eps) then
            if (y < brdy1 - eps) then
              freeflag = 1.
              if (size(bY1) < nbnewY1) then
                call resize(bY1,size(bY1),size(bY1)*2)
              end if
              bY1(nbnewY1) = n
              nbnewY1 = nbnewY1 + 1
            ! else if (y > brdy2 - iby * spy + eps) then
            else if (y > brdy2 + eps) then
              freeflag = 1.
              if (size(bY2) < nbnewY2) then
                call resize(bY2,size(bY2),size(bY2)*2)
              end if
              bY2(nbnewY2) = n
              nbnewY2 = nbnewY2 + 1
            end if
            if (dim == 3) then
              ! if (z < brdz1 + ibz * spz - eps) then
              if (z < brdz1 - eps) then
                freeflag = 1.
                if (size(bZ1) < nbnewZ1) then
                  call resize(bZ1,size(bZ1),size(bZ1)*2)
                end if
                bZ1(nbnewZ1) = n
                nbnewZ1 = nbnewZ1 + 1
              ! else if (z > brdz2 - ibz * spz + eps) then
              else if (z > brdz2 + eps) then
                freeflag = 1.
                if (size(bZ2) < nbnewZ2) then
                  call resize(bZ2,size(bZ2),size(bZ2)*2)
                end if
                bZ2(nbnewZ2) = n
                nbnewZ2 = nbnewZ2 + 1
              end if
            end if
          end if
          ptsz = size(pos, dim=2)
          if (ptsz < n) then
            call resize3r(pos, ptsz, ptsz*2)
            call resize(ptype, ptsz, ptsz*2)
          end if
          pos(1,n) = x
          pos(2,n) = y
          pos(3,n) = z
          ! if freeflag = 0 -> internal particle ()=> 0, else it is ghost ()=> 1
          ptype(n) = merge(1, 0, freeflag == 0)
          freenumber = merge(freenumber + 1, freenumber, freeflag == 0)
          n = n + 1
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
    call resize(ptype,n,n)
    write(*, "(A, F7.5, A, F7.5)") " # #   actual dx:   dx1=", pspc1, "  dx2=", pspc2
    write(*, "(A, I5, A, I5, A, I5)") " # #  dir.layers:   nx=", ix+1, &
             "   ny=", iy+1, "   nz=", iz+1
    write(*, "(A, I8, A, I8, A)") " # #    p.number:  ", n, &
              " total   ", freenumber, " real"
    print *, '# #    border-x:', brdx1, brdx2
    print *, '# #    border-y:', brdy1, brdy2
    print *, '# #    border-z:', brdz1, brdz2
    print *, '# #   № bd.pt X:', nbnewX1, nbnewX2
    print *, '# #   № bd.pt Y:', nbnewY1, nbnewY2
    print *, '# #   № bd.pt Z:', nbnewZ1, nbnewZ2

    call set_particles_numbers(n, abs(nb))
    call set_border(11, nbnewX1, bX1)
    call set_border(12, nbnewX2, bX2)
    call set_border(21, nbnewY1, bY1)
    call set_border(22, nbnewY2, bY2)
    call set_border(31, nbnewZ1, bZ1)
    call set_border(32, nbnewZ2, bZ2)

  end subroutine place_uniform

  subroutine place_close_packed_fcc(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, nb, pos)
    real, allocatable, intent(inout) :: pos(:,:)
    real, intent(inout)  :: pspc1, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2
    integer, allocatable :: bX1(:), bY1(:), bZ1(:), bX2(:), bY2(:), bZ2(:), bxprev(:), bxcur(:)
    integer              :: dim, nb, bdx, bdy, bdz, freeflag, freenumber, &
                            i, j, k, nx, ny, nz, n, nbnewX1, nbnewY1, nbnewZ1, nbnewX2, nbnewY2, nbnewZ2
    real                 :: dx, dy, dz, x, y, z, eps, sfxy, sfxz, sfyz, cy, sfbdx2
    eps = epsilon(0.)
    call get_dim(dim)

    allocate(bX1(1))
    allocate(bX2(0))
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

    nx = int((brdx2-brdx1)/pspc1)
    dx = merge(0., (brdx2-brdx1)/nx, nx == 0) ! Delta x
    dy = merge(0., sqrt(3./4.)*dx,   nx == 0) ! h
    cy = (brdy2+brdy1)/2
    ny = merge(0, int((brdy2 - cy)/dy), dy == 0)
    brdy2 = cy + ny * dy
    brdy1 = cy - ny * dy
    ny = 2*ny
    dz = merge(0., sqrt(2./3.)*dx,   nx == 0) ! delta z
    nz = merge(0, int((brdz2-brdz1)/dz), dz == 0)
    pspc1 = dx

    sfxy = 0.5   * dx ! shift on 'x' for the next y_row
    sfxz = 0.5   * dx ! shift on 'x' for the next z_row
    sfyz = 1./3. * dy ! shift on 'x' for the next z_row

    bdx = nb
    bdy = merge(nb, 0, dim > 1)
    bdz = merge(nb, 0, dim == 3)
    allocate(bxprev(ny+1+2*bdy))
    allocate(bxcur(ny+1+2*bdy))
    bxcur(:) = 0

    call set_sqare_box_sides(nx+1+2*bdx, ny+1+2*bdy, nz+1+2*bdz)
    freenumber = 0
    do i = (0-bdx),(nx+bdx)
      x = brdx1 + i * dx
      do j = (0-bdy),(ny+bdy)
        y = brdy1 + j * dy
        do k = (0-bdz),(nz+bdz)
          freeflag = 0
          if ((i /= (nx+bdx)).or.(mod(j,2)==0)) then
            sfbdx2 = 0.
            if (mod(j,2)==0) then
              sfbdx2 = 0.
            else
              sfbdx2 = dx
            end if
            z = brdz1 + k * dz
            if (x < brdx1 - eps) then
              freeflag = 1.
              if (size(bX1) < nbnewX1) then
                call resize(bX1,size(bX1),size(bX1)*2)
              end if
              bX1(nbnewX1) = n
              nbnewX1 = nbnewX1 + 1
            else if (x > brdx2 - sfbdx2 + eps) then
              freeflag = 1.
              ! if (size(bX2) < nbnewX2) then
              !   call resize(bX2,size(bX2),size(bX2)*2)
              ! end if
              ! bX2(nbnewX2) = n
              bxcur(nbnewX2) = n
              nbnewX2 = nbnewX2 + 1
            end if
            if (dim > 1) then
              if (y < brdy1 - eps) then
                ! print *, "BRDY1"
                freeflag = 1.
                if (size(bY1) < nbnewY1) then
                  call resize(bY1,size(bY1),size(bY1)*2)
                end if
                bY1(nbnewY1) = n
                nbnewY1 = nbnewY1 + 1
              else if (y > brdy2 + eps) then
                ! print *, "BRDY2"
                freeflag = 1.
                if (size(bY2) < nbnewY2) then
                  call resize(bY2,size(bY2),size(bY2)*2)
                end if
                bY2(nbnewY2) = n
                nbnewY2 = nbnewY2 + 1
              end if
              if (dim == 3) then
                if (z < brdz1 - eps) then
                  freeflag = 1.
                  if (size(bZ1) < nbnewZ1) then
                    call resize(bZ1,size(bZ1),size(bZ1)*2)
                  end if
                  bZ1(nbnewZ1) = n
                  nbnewZ1 = nbnewZ1 + 1
                else if (z > brdz2 + eps) then
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
            ! to set borders like for uniform and then bend them
            if (mod(j,2)==0) then
              pos(1,n) = x
            else
              pos(1,n) = x+sfxy
            end if
            pos(2,n) = y
            pos(3,n) = z
            n = n + 1
            freenumber = merge(freenumber + 1, freenumber, freeflag == 0)
          end if
        end do
      end do
      if (i >= nx) then
        if (i == nx) then
          bxprev(:) = bxcur(:)
          nbnewX2 = 1
          bxcur(:) = 0
        else if (i == nx+bdx) then
          nbnewX2 = size(bX2)
          call resize(bX2,nbnewX2,nbnewX2+size(bxcur))
          do k = 1,size(bxcur)
            if (mod(k,2)==1) then
              bX2(nbnewX2+k) = bxcur(k/2+1)
            else
              bX2(nbnewX2+k) = bxprev(k/2)
            end if
          end do
          nbnewX2 = size(bX2) + 1
          bxcur(:) = 0
          deallocate(bxcur)
          deallocate(bxprev)
        else
          ! In fact old size of border particles on x 2
          nbnewX2 = size(bX2)
          call resize(bX2,nbnewX2,nbnewX2+size(bxcur))
          do k = 1,size(bxcur)
            if (mod(k,2)==1) then
              bX2(nbnewX2+k) = bxcur(k)
            else
              bX2(nbnewX2+k) = bxprev(k/2)
              bxprev(k/2) = bxcur(k)
            end if
          end do
          nbnewX2 = 1
          bxcur(:) = 0
        end if
      end if
    end do

    nbnewX1 = nbnewX1 - 1
    nbnewY1 = nbnewY1 - 1
    nbnewZ1 = nbnewZ1 - 1
    nbnewX2 = nbnewX2 - 1
    nbnewY2 = nbnewY2 - 1
    nbnewZ2 = nbnewZ2 - 1
    n = n - 1

    call resize3r(pos,n,n)

    write(*, "(A, F7.5, A, F7.5, A, F7.5)") " # # hex.spacing:   dx=", pspc1, " hy=", dy, " hz=", dz
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
  end subroutine place_close_packed_fcc
end module initpositions
