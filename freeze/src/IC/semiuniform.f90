module semiuniform
  use kernel
  use utils, only:resize
  use BC

  implicit none

  public :: make_semiuniform

  private

contains

  subroutine make_semiuniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos)
    real, allocatable, intent(inout) :: pos(:,:)
    real, intent(in)     :: brdx1, brdx2, brdy1, brdy2, brdz1, brdz2
    real, intent(inout)  :: pspc1, pspc2
    integer, allocatable :: bX1(:), bY1(:), bZ1(:), bX2(:), bY2(:), bZ2(:)
    integer              :: i, dim, nb, freeflag, freenumber, &
                            n, nbnewX1, nbnewY1, nbnewZ1, nbnewX2, nbnewY2, nbnewZ2
    real                 :: x, y, z, sp, eps

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

    freenumber = 0

    eps = 1e-5

    x = brdx1! - pspc1 * nbx
    ! do while ((x >= brdx1 - pspc1 * nbx).and.(x <= brdx2 + pspc2 * nbx))
    do while ((x >= brdx1).and.(x <= brdx2 + eps))
      if (x < 0) then
        sp = pspc1
      else
        sp = pspc2
        if ((x + sp > brdx2) .and. (brdx2 - x > pspc2/2)) then
          sp = brdx2 - x
        end if
      end if
      y = brdy1! - pspc1 * nby
      ! do while ((y >= brdy1 - pspc1 * nby).and.(y <= brdy2 + pspc2 * nby))
      do while ((y >= brdy1).and.(y <= brdy2 + eps))
        z = brdz1! - pspc1 * nbz
        ! do while ((z >= brdz1 - pspc1 * nbz).and.(z <= brdz2 + pspc2 * nbz))
        do while ((z >= brdz1).and.(z <= brdz2 + eps))
          freeflag = 0
          if (x < brdx1 + nb * sp) then
            freeflag = 1.
            if (size(bX1) < nbnewX1) then
              call resize(bX1,size(bX1),size(bX1)*2)
            end if
            bX1(nbnewX1) = n
            nbnewX1 = nbnewX1 + 1
          else if (x > brdx2 - nb * sp - eps) then
            freeflag = 1.
            if (size(bX2) < nbnewX2) then
              call resize(bX2,size(bX2),size(bX2)*2)
            end if
            bX2(nbnewX2) = n
            nbnewX2 = nbnewX2 + 1
          end if
          if (dim > 1) then
            ! if (y < brdy1) then
            if (y < brdy1 + nb * sp) then
              freeflag = 1.
              if (size(bY1) < nbnewY1) then
                call resize(bY1,size(bY1),size(bY1)*2)
              end if
              bY1(nbnewY1) = n
              nbnewY1 = nbnewY1 + 1
            ! else if (y > brdy2) then
          else if (y > brdy2 - nb * sp - eps) then
              freeflag = 1.
              if (size(bY2) < nbnewY2) then
                call resize(bY2,size(bY2),size(bY2)*2)
              end if
              bY2(nbnewY2) = n
              nbnewY2 = nbnewY2 + 1
            end if
            if (dim == 3) then
              if (z < brdz1 + nb * sp) then
                freeflag = 1.
                if (size(bZ1) < nbnewZ1) then
                  call resize(bZ1,size(bZ1),size(bZ1)*2)
                end if
                bZ1(nbnewZ1) = n
                nbnewZ1 = nbnewZ1 + 1
              else if (z > brdz2 - nb * sp - eps) then
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
            call resize(pos,size(pos, dim=2),size(pos, dim=2)*2)
          end if
          pos(1,n) = x
          pos(2,n) = y
          pos(3,n) = z
          n = n + 1
          freenumber = merge(freenumber + 1, freenumber, freeflag == 0)
          z = z + sp
        end do
        y = y + sp
      end do
      x = x + sp
    end do

    do i=1,n-1
    end do

    nbnewX1 = nbnewX1 - 1
    nbnewY1 = nbnewY1 - 1
    nbnewZ1 = nbnewZ1 - 1
    nbnewX2 = nbnewX2 - 1
    nbnewY2 = nbnewY2 - 1
    nbnewZ2 = nbnewZ2 - 1
    n = n - 1

    call resize(pos,n,n)

    write(*, "(A, F7.5, A, F7.5)") " # #      actual dx:   x1=", pspc1, "   x2=", pspc2
    print *, '# #         placed:', n
    print *, '# #     freenumber:', freenumber
    print *, '# #       border-x:', nbnewX1, nbnewX2
    print *, '# #       border-y:', nbnewY1, nbnewY2
    print *, '# #       border-z:', nbnewZ1, nbnewZ2

    call set_particles_numbers(n, abs(nb))
    call set_border(11, nbnewX1, bX1)
    call set_border(12, nbnewX2, bX2)
    call set_border(21, nbnewY1, bY1)
    call set_border(22, nbnewY2, bY2)
    call set_border(31, nbnewZ1, bZ1)
    call set_border(32, nbnewZ2, bZ2)

  end subroutine make_semiuniform
end module semiuniform
