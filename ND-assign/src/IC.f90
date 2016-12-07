module IC
  use kernel
  use BC

  implicit none

  public :: setup

  private
  real, parameter     :: pi = 4.*atan(1.)

contains

  subroutine setup(tt, kt, dim, n, sk, g, cv, pspc1, pspc2, &
    pos, vel, acc, mas, den, sln, prs, iu, du, cf, kcf, dcf)
    character(len=*), intent(in) :: tt, kt
    integer, intent(in)  :: dim
    real, intent(in)     :: sk, g, cv
    real, allocatable, intent(inout) :: pos(:,:), vel(:,:), acc(:,:), &
                                        mas(:), den(:), sln(:), prs(:), iu(:), du(:), cf(:), kcf(:), dcf(:)
    real, intent(inout)  :: pspc1, pspc2
    integer, intent(out) :: n
    real                 :: kr, prs1, prs2, rho1, rho2, k1, k2, x, y, z, sp, spx, spy, spz, &
                            brdx1, brdx2, brdy1, brdy2, brdz1, brdz2
    integer              :: i, j, k, nb, ix, iy, iz, bdx, bdy, bdz, ibx, iby, ibz, &
                            nbnewX1, nbnewY1, nbnewZ1, nbnewX2, nbnewY2, nbnewZ2, freeflag, freenumber
    integer, allocatable :: bX1(:), bY1(:), bZ1(:), bX2(:), bY2(:), bZ2(:)

    allocate(bX1(1))
    allocate(bX2(1))
    allocate(bY1(1))
    allocate(bY2(1))
    allocate(bZ1(1))
    allocate(bZ2(1))
    allocate(pos(3,1))

    call set_dim(dim)
    call set_tasktype(tt)
    call set_kerntype(kt)
    call get_krad(kr)

    nb = 0
    brdx1 = 0
    brdx2 = 0
    brdy1 = 0
    brdy2 = 0
    brdz1 = 0
    brdz2 = 0
    rho1 = 1.
    rho2 = 1.
    prs1 = 1.
    prs2 = 1.
    k1 = 1.
    k2 = 1.

    select case (tt)
    case ('hydroshock')
      nb = 4
      prs1 = 1.
      prs2 = 0.1
      rho1 = 1.
      rho2 = 0.125
    case ('infslb')
      nb = 1
    case ('hc-sinx')
      nb = int(kr * sk) + 1
      brdx1 = -1.
      brdx2 = 1.
      if (dim > 1) then
        brdy1 = -0.25
        brdy2 = 0.25
      else
        brdy1 = 0.
        brdy2 = 0.
      end if
      if (dim == 3) then
        brdz1 = -0.25
        brdz2 = 0.25
      else
        brdz1 = 0.
        brdz2 = 0.
      end if
    end select

    n = 1
    nbnewX1 = 1
    nbnewY1 = 1
    nbnewZ1 = 1
    nbnewX2 = 1
    nbnewY2 = 1
    nbnewZ2 = 1

    if (pspc1 /= pspc2) then
      ! ! if (dim > 1) then
      ! !   nby = nb
      ! !   if (dim == 3) then
      ! !     nbz = nb
      ! !   end if
      ! ! end if
      ! x = brdx1! - pspc1 * nbx
      ! ! do while ((x >= brdx1 - pspc1 * nbx).and.(x <= brdx2 + pspc2 * nbx))
      ! do while ((x >= brdx1).and.(x <= brdx2))
      !   if (x < 0) then
      !     sp = pspc1
      !   else
      !     sp = pspc2
      !   end if
      !   y = brdy1! - pspc1 * nby
      !   ! do while ((y >= brdy1 - pspc1 * nby).and.(y <= brdy2 + pspc2 * nby))
      !   do while ((y >= brdy1).and.(y <= brdy2))
      !     z = brdz1! - pspc1 * nbz
      !     ! do while ((z >= brdz1 - pspc1 * nbz).and.(z <= brdz2 + pspc2 * nbz))
      !     do while ((z >= brdz1).and.(z <= brdz2))
      !       pos(1,n) = x
      !       pos(2,n) = y
      !       pos(3,n) = z
      !       if (x < brdx1 + nb * sp) then
      !         bX1(nbnewX1) = n
      !         nbnewX1 = nbnewX1 + 1
      !       else if (x > brdx2 - nb * sp) then
      !         bX2(nbnewX2) = n
      !         nbnewX2 = nbnewX2 + 1
      !       end if
      !       if (dim > 1) then
      !         ! if (y < brdy1) then
      !         if (y < brdy1 + nb * sp) then
      !           bY1(nbnewY1) = n
      !           nbnewY1 = nbnewY1 + 1
      !         ! else if (y > brdy2) then
      !         else if (y > brdy2 - nb * sp) then
      !           bY2(nbnewY2) = n
      !           nbnewY2 = nbnewY2 + 1
      !         end if
      !         if (dim == 3) then
      !           if (z < brdz1 + nb * sp) then
      !             bZ1(nbnewZ1) = n
      !             nbnewZ1 = nbnewZ1 + 1
      !           else if (z > brdz2 - nb * sp) then
      !             bZ2(nbnewZ2) = n
      !             nbnewZ2 = nbnewZ2 + 1
      !           end if
      !         end if
      !       end if
      !       z = z + sp
      !       n = n + 1
      !     end do
      !     y = y + sp
      !   end do
      !   x = x + sp
      ! end do
    else
      ix = int((brdx2-brdx1)/pspc1)
      spx = merge(0.,(brdx2-brdx1)/ix, ix == 0)
      if (dim > 1) then
        brdy1 = - int(kr * sk) * spx * 2
        brdy2 =   int(kr * sk) * spx * 2
      else
        brdy1 = 0.
        brdy2 = 0.
      end if
      if (dim == 3) then
        brdz1 = - int(kr * sk) * spx * 2
        brdz2 =   int(kr * sk) * spx * 2
      else
        brdz1 = 0.
        brdz2 = 0.
      end if
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
    end if

    nbnewX1 = nbnewX1 - 1
    nbnewY1 = nbnewY1 - 1
    nbnewZ1 = nbnewZ1 - 1
    nbnewX2 = nbnewX2 - 1
    nbnewY2 = nbnewY2 - 1
    nbnewZ2 = nbnewZ2 - 1
    n = n - 1

    call resize3r(pos,n,n)
    allocate(vel(3,n))
    allocate(acc(3,n))
    allocate(mas(n))
    allocate(den(n))
    allocate(sln(n))
    allocate(prs(n))
    allocate(iu(n))
    allocate(du(n))
    allocate(cf(n))
    allocate(kcf(n))
    allocate(dcf(n))

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
    !
    ! common values
    !
    do i=1,n
      sp = merge(pspc1, pspc2, pos(1,i) < 0)

      vel(:,i) = 0.
      acc(:,i) = 0.
      dcf(i) = 0.
      sln(i) = sk * sp
      prs(i) = 0

      if (pos(1,i) < 0) then
        mas(i) = (sp**dim) * rho1
        den(i) = rho1
        prs(i) = prs1
        kcf(i) = k1
      else
        mas(i) = (sp**dim) * rho2
        den(i) = rho2
        prs(i) = prs2
        kcf(i) = k2
      end if

      select case (tt)
      case ('hydroshock')
        iu(i) = merge(prs1/(g-1)/rho1, prs2/(g-1)/rho2, pos(1,i) < 0)
      case ('infslb')
        cf(i) = merge(0., 1., pos(1,i) < 0)
        iu(i) = cf(i) / cv
      case ('hc-sinx')
        cf(i)  = sin(2 * pi * (pos(1,i) + 1.) / abs(brdx2-brdx1))
        iu(i) = cf(i) / cv
      end select
    end do
  end subroutine setup

  subroutine resize(array, oldsize, newsize)
    integer, intent(in)                 :: newsize, oldsize
    integer, intent(inout), allocatable :: array(:)
    integer, allocatable                :: tmp(:)
    integer                             :: i

    allocate(tmp(newsize))
    do i=1,oldsize
      tmp(i) = array(i)
    end do
    deallocate(array)
    allocate(array(newsize))
    do i=1,oldsize
      array(i) = tmp(i)
    end do
  end subroutine

  subroutine resize3r(array, oldsize, newsize)
    integer, intent(in)              :: newsize, oldsize
    real, intent(inout), allocatable :: array(:,:)
    real, allocatable                :: tmp(:,:)
    integer                          :: i

    allocate(tmp(3,newsize))
    do i=1,oldsize
      tmp(:,i) = array(:,i)
    end do
    deallocate(array)
    allocate(array(3,newsize))
    do i=1,oldsize
      array(:,i) = tmp(:,i)
    end do
  end subroutine
end module IC
