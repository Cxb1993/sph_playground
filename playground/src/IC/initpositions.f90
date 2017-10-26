module InitPositions
  use const
  use kernel
  use state,  only: getdim
  use ArrayResize,  only: resize
  use BC

  implicit none

  public :: uniform, semiuniform, place_close_packed_fcc

  private

contains

  subroutine uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos, ptype)
    real,    allocatable, intent(inout) :: pos(:,:)
    integer, allocatable, intent(inout) :: ptype(:)
    real, intent(inout)  :: pspc1, pspc2, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2

    integer, allocatable :: bX1(:), bY1(:), bZ1(:), bX2(:), bY2(:), bZ2(:)
    integer              :: dim, nb, bdx, bdy, bdz, isborder, freenumber, &
                            i, j, k, ix, iy, iz, n, nbnewX1, nbnewY1, nbnewZ1, nbnewX2, nbnewY2, nbnewZ2,&
                            ptsz !particle array size
                            ! , ibx, iby, ibz ! number thet were used for 'niternal' border setting
    real                 :: spx, spy, spz, x, y, z, eps

    ! eps = 10*epsilon(0.)
    eps = eps0
    eps = -eps0

    call getdim(dim)

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
          isborder = 0
          z = brdz1 + k * spz
          ! if (x < brdx1 + ibx * spx - eps) then
          if (x < brdx1 + eps) then
          ! if (x < brdx1 - eps) then
            isborder = 1.
            if (size(bX1) < nbnewX1) then
              call resize(bX1,size(bX1),size(bX1)*2)
            end if
            bX1(nbnewX1) = n
            nbnewX1 = nbnewX1 + 1
          ! else if (x > brdx2 - ibx * spx + eps) then
          else if (x > brdx2 - eps) then
          ! else if (x > brdx2 + eps) then
            isborder = 1.
            if (size(bX2) < nbnewX2) then
              call resize(bX2,size(bX2),size(bX2)*2)
            end if
            bX2(nbnewX2) = n
            nbnewX2 = nbnewX2 + 1
          end if
          if (dim > 1) then
            ! if (y < brdy1 + iby * spy - eps) then
            if (y < brdy1 + eps) then
            ! if (y < brdy1 - eps) then
              isborder = 1.
              if (size(bY1) < nbnewY1) then
                call resize(bY1,size(bY1),size(bY1)*2)
              end if
              bY1(nbnewY1) = n
              nbnewY1 = nbnewY1 + 1
            ! else if (y > brdy2 - iby * spy + eps) then
            else if (y > brdy2 - eps) then
            ! else if (y > brdy2 + eps) then
              isborder = 1.
              if (size(bY2) < nbnewY2) then
                call resize(bY2,size(bY2),size(bY2)*2)
              end if
              bY2(nbnewY2) = n
              nbnewY2 = nbnewY2 + 1
            end if
            if (dim == 3) then
              ! if (z < brdz1 + ibz * spz - eps) then
              if (z < brdz1 + eps) then
              ! if (z < brdz1 - eps) then
                isborder = 1.
                if (size(bZ1) < nbnewZ1) then
                  call resize(bZ1,size(bZ1),size(bZ1)*2)
                end if
                bZ1(nbnewZ1) = n
                nbnewZ1 = nbnewZ1 + 1
              ! else if (z > brdz2 - ibz * spz + eps) then
              else if (z > brdz2 - eps) then
              ! else if (z > brdz2 + eps) then
                isborder = 1.
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
            call resize(pos, ptsz, ptsz*2)
            call resize(ptype, ptsz, ptsz*2)
          end if
          pos(1,n) = x
          pos(2,n) = y
          pos(3,n) = z
          ! if isborder = 0 -> internal particle ()=> 0, else it is ghost ()=> 1
          ptype(n) = merge(1, 0, isborder == 0)
          freenumber = merge(freenumber + 1, freenumber, isborder == 0)
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

    call resize(pos,n,n)
    call resize(ptype,n,n)
    write(*, "(A, F7.5, A, F7.5)") " # #        actual dx:   dx1=", pspc1, "  dx2=", pspc2
    write(*, "(A, I5, A, I5, A, I5)") " # #       dir.layers:   nx=", ix+1, &
             "   ny=", iy+1, "   nz=", iz+1
    write(*, "(A, I16, A, I16, A)") " # #         p.number:  ", n, &
              " total   ", freenumber, " real"
    print *, '# #         border-x:', brdx1, brdx2
    print *, '# #         border-y:', brdy1, brdy2
    print *, '# #         border-z:', brdz1, brdz2
    print *, '# #        № bd.pt X:', nbnewX1, nbnewX2
    print *, '# #        № bd.pt Y:', nbnewY1, nbnewY2
    print *, '# #        № bd.pt Z:', nbnewZ1, nbnewZ2

    call set_particles_numbers(n, abs(nb))
    call resize(bx1,nbnewX1,nbnewX1)
    call resize(bx2,nbnewX2,nbnewX2)
    call resize(by1,nbnewY1,nbnewY1)
    call resize(by2,nbnewY2,nbnewY2)
    call resize(bz1,nbnewZ1,nbnewZ1)
    call resize(bz2,nbnewZ2,nbnewZ2)
    call setBorder(11, bX1)
    call setBorder(12, bX2)
    call setBorder(21, bY1)
    call setBorder(22, bY2)
    call setBorder(31, bZ1)
    call setBorder(32, bZ2)
  end subroutine

  subroutine semiuniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos, ptype)
    use list

    real, intent(in)    :: brdx1, brdx2, brdy1, brdy2, brdz1, brdz2

    real, allocatable, intent(inout)    :: pos(:,:)
    integer, allocatable, intent(inout) :: ptype(:)
    real, intent(inout) :: pspc1, pspc2

    type(intlist)       :: bx1, bx2, by1, by2, bz1, bz2
    type(intlist)       :: inbx1, inbx2, inby1, inby2, inbz1, inbz2
    integer             :: dim, nb, isborder, freenumber, n, ptsz, &
                            ix1, ix2, iy1, iz1, bx, by, bz, i, j, k
    real                :: sp, eps

    ! eps = -eps0
    eps = eps0

    call getdim(dim)

    allocate(pos(3,1))
    allocate(ptype(1))

    n = 1

    freenumber = 0

    ix1 = int(-brdx1/pspc1)
    ix2 = int(brdx2/pspc2)
    ! dp1 = merge(0.,(-brdx1)/ix1, ix1 == 0)
    ! dp2 = merge(0.,(brdx2)/ix2, ix2 == 0)
    ! pspc1 = dp1
    ! pspc2 = dp2
    bx = nb
    by = merge(0, nb, abs(brdy1-brdy2) < eps)
    bz = merge(0, nb, abs(brdz1-brdz2) < eps)
    ! print*, (-bx -ix1), (ix2 +bx), dp1*(-bx -ix1), dp2*(ix2 +bx)
    ! read*
    do i = (-bx -ix1), (ix2 +bx)
      if (i <= 0) then
        sp = pspc1
      else
        sp = pspc2
      end if
      iy1 = int((brdy2-brdy1)/sp/2)
      iz1 = int((brdz2-brdz1)/sp/2)
      ! print*,(-by -iy1), (iy1 +by), by, iy1
      ! print*, '------'
      do j = (-by -iy1), (iy1 +by)
        ! print*, j
        do k = (-bz -iz1), (iz1 +bz)
          isborder = 0
          if (i < -ix1) then
            isborder = 1
            call bx1%append(n)
          else if (i < -ix1 + bx) then
            call inbx1%append(n)
          else if (i <= ix2 - bx) then
          else if (i <= ix2) then
            call inbx2%append(n)
          else
            isborder = 1
            call bx2%append(n)
          end if
          if (dim > 1) then
            if (j < -iy1) then
              isborder = 1
              call by1%append(n)
            else if (j < -iy1 + by) then
              call inby1%append(n)
            else if (j <= iy1 - by) then
            else if (j <= iy1) then
              call inby2%append(n)
            else
              isborder = 1
              call by2%append(n)
            end if
            if (dim == 3) then
              if (k < -iz1) then
                isborder = 1
                call bz1%append(n)
              else if (k < -iz1 + bz) then
                call inbz1%append(n)
              else if (k <= iz1 - bz) then
              else if (k <= iz1) then
                call inbz2%append(n)
              else
                isborder = 1
                call bz2%append(n)
              end if
            end if
          end if
          ptsz = size(pos, dim=2)
          if (ptsz < n) then
            call resize(pos, ptsz, ptsz*2)
            call resize(ptype, ptsz, ptsz*2)
          end if
          pos(1,n) = sp*i
          pos(2,n) = sp*j
          pos(3,n) = sp*k
          freenumber = merge(freenumber + 1, freenumber, isborder == 0)
          ptype(n) = merge(1, 0, isborder == 0)
          n = n + 1
        end do
      end do
      ! read*
    end do

    n = n - 1

    call resize(pos,n,n)
    call resize(ptype,n,n)

    ! print*, by1%toarr()
    ! print*, '--------------------'
    ! print*, by2%toarr()
    ! print*, '--------------------'
    ! print*, inby1%toarr()
    ! print*, '--------------------'
    ! print*, inby2%toarr()
    ! print*, '--------------------'
    ! read*

    write(*, "(A, F7.5, A, F7.5)") " # #        actual dx:   x1=", pspc1, "   x2=", pspc2
    write(*, "(A, I16, A, I16, A)") " # #         p.number:  ", n, &
              " total   ", freenumber, " real"
    print *, '# #       border x:', bx1%llen(), bx2%llen()
    print *, '# # inner border x:', inbx1%llen(), inbx2%llen()
    print *, '# #       border y:', by1%llen(), by2%llen()
    print *, '# # inner border y:', inby1%llen(), inby2%llen()
    print *, '# #       border z:', bz1%llen(), bz2%llen()
    print *, '# # inner border z:', inbz1%llen(), inbz2%llen()

    call set_particles_numbers(n, abs(nb))

    call setBorder(11, bx1%toarr())
    call setBorder(12, bx2%toarr())
    call setBorder(21, by1%toarr())
    call setBorder(22, by2%toarr())
    call setBorder(31, bz1%toarr())
    call setBorder(32, bz2%toarr())

    call setBorderInside(11, inbx1%toarr())
    call setBorderInside(12, inbx2%toarr())
    call setBorderInside(21, inby1%toarr())
    call setBorderInside(22, inby2%toarr())
    call setBorderInside(31, inbz1%toarr())
    call setBorderInside(32, inbz2%toarr())

    ! print*, bY1
    ! read*
  end subroutine
  subroutine place_close_packed_fcc(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, nb, pos)
    real, allocatable, intent(inout) :: pos(:,:)
    real, intent(inout)  :: pspc1, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2
    integer, allocatable :: bX1(:), bY1(:), bZ1(:), bX2(:), bY2(:), bZ2(:), bxprev(:), bxcur(:)
    integer              :: dim, nb, bdx, bdy, bdz, isborder, freenumber, &
                            i, j, k, nx, ny, nz, n, nbnewX1, nbnewY1, nbnewZ1, nbnewX2, nbnewY2, nbnewZ2
    real                 :: dx, dy, dz, x, y, z, eps, sfxy, sfxz, sfyz, cy, sfbdx2
    eps = epsilon(0.)
    call getdim(dim)

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
    ny = merge(0, int((brdy2 - cy)/dy), nx == 0)
    brdy2 = cy + ny * dy
    brdy1 = cy - ny * dy
    ny = 2*ny
    dz = merge(0., sqrt(2./3.)*dx,   nx == 0) ! delta z
    nz = merge(0, int((brdz2-brdz1)/dz), nx == 0)
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
          isborder = 0
          if ((i /= (nx+bdx)).or.(mod(j,2)==0)) then
            sfbdx2 = 0.
            if (mod(j,2)==0) then
              sfbdx2 = 0.
            else
              sfbdx2 = dx
            end if
            z = brdz1 + k * dz
            if (x < brdx1 - eps) then
              isborder = 1.
              if (size(bX1) < nbnewX1) then
                call resize(bX1,size(bX1),size(bX1)*2)
              end if
              bX1(nbnewX1) = n
              nbnewX1 = nbnewX1 + 1
            else if (x > brdx2 - sfbdx2 + eps) then
              isborder = 1.
              ! if (size(bX2) < nbnewX2) then
              !   call iresize(bX2,size(bX2),size(bX2)*2)
              ! end if
              ! bX2(nbnewX2) = n
              bxcur(nbnewX2) = n
              nbnewX2 = nbnewX2 + 1
            end if
            if (dim > 1) then
              if (y < brdy1 - eps) then
                ! print *, "BRDY1"
                isborder = 1.
                if (size(bY1) < nbnewY1) then
                  call resize(bY1,size(bY1),size(bY1)*2)
                end if
                bY1(nbnewY1) = n
                nbnewY1 = nbnewY1 + 1
              else if (y > brdy2 + eps) then
                ! print *, "BRDY2"
                isborder = 1.
                if (size(bY2) < nbnewY2) then
                  call resize(bY2,size(bY2),size(bY2)*2)
                end if
                bY2(nbnewY2) = n
                nbnewY2 = nbnewY2 + 1
              end if
              if (dim == 3) then
                if (z < brdz1 - eps) then
                  isborder = 1.
                  if (size(bZ1) < nbnewZ1) then
                    call resize(bZ1,size(bZ1),size(bZ1)*2)
                  end if
                  bZ1(nbnewZ1) = n
                  nbnewZ1 = nbnewZ1 + 1
                else if (z > brdz2 + eps) then
                  isborder = 1.
                  if (size(bZ2) < nbnewZ2) then
                    call resize(bZ2,size(bZ2),size(bZ2)*2)
                  end if
                  bZ2(nbnewZ2) = n
                  nbnewZ2 = nbnewZ2 + 1
                end if
              end if
            end if
            if (size(pos, dim=2) < n) then
              call resize(pos,size(pos, dim=2),size(pos, dim=2)*2)
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
            freenumber = merge(freenumber + 1, freenumber, isborder == 0)
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

    call resize(pos,n,n)

    write(*, "(A, F7.5, A, F7.5, A, F7.5)") " # # hex.spacing:   dx=", pspc1, " hy=", dy, " hz=", dz
    print *, '# #      placed:', n
    print *, '# #  freenumber:', freenumber
    print *, '# #    border-x:', nbnewX1, nbnewX2
    print *, '# #    border-y:', nbnewY1, nbnewY2
    print *, '# #    border-z:', nbnewZ1, nbnewZ2

    call set_particles_numbers(n, abs(nb))
    call resize(bx1,nbnewX1,nbnewX1)
    call resize(bx2,nbnewY1,nbnewY1)
    call resize(by1,nbnewZ1,nbnewZ1)
    call resize(by2,nbnewX2,nbnewX2)
    call resize(bz1,nbnewY2,nbnewY2)
    call resize(bz2,nbnewZ2,nbnewZ2)
    call setBorder(11, bX1)
    call setBorder(12, bX2)
    call setBorder(21, bY1)
    call setBorder(22, bY2)
    call setBorder(31, bZ1)
    call setBorder(32, bZ2)
  end subroutine
end module
