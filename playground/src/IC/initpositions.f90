module InitPositions
  use const
  use kernel
  use state,  only: getdim
  use ArrayResize,  only: resize
  use BC

  implicit none

  public :: uniform

  private

contains
  subroutine uniform(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos, ptype, randomise)
    use list
    use NeighbourSearch, only: FindPeriodicBorder

    real, intent(in)    :: brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, randomise
    real, allocatable, intent(inout)    :: pos(:,:)
    integer, allocatable, intent(inout) :: ptype(:)
    real, intent(inout)  :: pspc1, pspc2

    type(intlist)        :: bx1, bx2, by1, by2, bz1, bz2
    integer, allocatable :: inbx1(:), inbx2(:), inby1(:), inby2(:), inbz1(:), inbz2(:)
    integer              :: dim, nb, isborder, freenumber, n, ptsz, &
                             ix1, ix2, iy1, iy2, iz1, iz2, bx, by, bz, i, j, k
    real                 :: sp, eps

    ! eps = -eps0
    eps = eps0

    call getdim(dim)

    allocate(pos(3,1))
    allocate(ptype(1))

    n = 1

    freenumber = 0

    print *, '# #       x in [',brdx1,":",brdx2,"]"
    print *, '# #       y in [',brdy1,":",brdy2,"]"
    print *, '# #       z in [',brdz1,":",brdz2,"]"

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
      iy1 = int(-brdy1/sp)
      iy2 = int(brdy2/sp)
      iz1 = int(-brdz1/sp)
      iz2 = int(brdz2/sp)
      ! print*, iy1, iy2, sp, iy1*sp, (iy1+1)*sp,
      ! read*
      ! print*,(-by -iy1), (iy1 +by), by, iy1
      ! print*, '------'
      do j = (-by -iy1), (iy2 +by)
        ! print*, j
        do k = (-bz -iz1), (iz2 +bz)
          isborder = 0
          if (i < -ix1) then
            isborder = 1
            call bx1%append(n)
          else if (i > ix2) then
            isborder = 1
            call bx2%append(n)
          end if
          if (dim > 1) then
            ! print*, '--------'
            ! print*, n, i, j, k
            if (j < -iy1) then
              isborder = 1
              call by1%append(n)
            else if (j > iy2) then
              isborder = 1
              call by2%append(n)
            end if
            if (dim == 3) then
              if (k < -iz1) then
                isborder = 1
                call bz1%append(n)
              else if (k > iz2) then
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
          pos(1,n) = sp*i + randomise*sp*(2*(rand() - 0.5))
          pos(2,n) = sp*j + randomise*sp*(2*(rand() - 0.5))
          pos(3,n) = sp*k + randomise*sp*(2*(rand() - 0.5))
          freenumber = merge(freenumber + 1, freenumber, isborder == 0)
          ptype(n) = merge(1, 0, isborder == 0)
          ! print*, pos(:,n), ptype(n)
          n = n + 1
        end do
      end do
      ! read*
    end do

    n = n - 1

    call resize(pos,n,n)
    call resize(ptype,n,n)

    call FindPeriodicBorder(pos, pspc1, brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, &
                            bx1, bx2, by1, by2, bz1, bz2, &
                            inbx1, inbx2, inby1, inby2, inbz1, inbz2)

    write(*, "(A, F7.5, A, F7.5)") " # #        actual dx:   x1=", pspc1, "   x2=", pspc2
    write(*, "(A, I16, A, I16, A)") " # #         p.number:  ", n, &
              " total   ", freenumber, " real"
    print *, '# #       border x:', bx1%llen(), bx2%llen()
    print *, '# # inner border x:', size(inbx1), size(inbx2)
    print *, '# #       border y:', by1%llen(), by2%llen()
    print *, '# # inner border y:', size(inby1), size(inby2)
    print *, '# #       border z:', bz1%llen(), bz2%llen()
    print *, '# # inner border z:', size(inbz1), size(inbz2)

    call set_particles_numbers(n, abs(nb))

    call setBorder(11, bx1%toarr())
    call setBorder(12, bx2%toarr())
    call setBorder(21, by1%toarr())
    call setBorder(22, by2%toarr())
    call setBorder(31, bz1%toarr())
    call setBorder(32, bz2%toarr())

    call setBorderInside(11, inbx1)
    call setBorderInside(12, inbx2)
    call setBorderInside(21, inby1)
    call setBorderInside(22, inby2)
    call setBorderInside(31, inbz1)
    call setBorderInside(32, inbz2)

  end subroutine

  ! subroutine place_close_packed_fcc(brdx1, brdx2, brdy1, brdy2, brdz1, brdz2, pspc1, pspc2, nb, pos, ptype)
  !   use list
  !
  !   real, intent(in)    :: brdx1, brdx2, brdy1, brdy2, brdz1, brdz2
  !
  !   real, allocatable, intent(inout)    :: pos(:,:)
  !   integer, allocatable, intent(inout) :: ptype(:)
  !   real, intent(inout) :: pspc1, pspc2
  !
  !   type(intlist)       :: bx1, bx2, by1, by2, bz1, bz2
  !   type(intlist)       :: inbx1, inbx2, inby1, inby2, inbz1, inbz2
  !   integer             :: dim, nb, isborder, freenumber, n, ptsz, &
  !                           ix1, ix2, iy1, iz1, bx, by, bz, i, j, k
  !   real                :: sp, eps
  !
  !   ! eps = -eps0
  !   eps = eps0
  !
  !   call getdim(dim)
  !
  !   allocate(pos(3,1))
  !   allocate(ptype(1))
  !
  !   n = 1
  !
  !   freenumber = 0
  !
  !   ix1 = int(-brdx1/pspc1)
  !   ix2 = int(brdx2/pspc2)
  !   ! dp1 = merge(0.,(-brdx1)/ix1, ix1 == 0)
  !   ! dp2 = merge(0.,(brdx2)/ix2, ix2 == 0)
  !   ! pspc1 = dp1
  !   ! pspc2 = dp2
  !   bx = nb
  !   by = merge(0, nb, abs(brdy1-brdy2) < eps)
  !   bz = merge(0, nb, abs(brdz1-brdz2) < eps)
  !   ! print*, (-bx -ix1), (ix2 +bx), dp1*(-bx -ix1), dp2*(ix2 +bx)
  !   ! read*
  !   do i = (-bx -ix1), (ix2 +bx)
  !     if (i <= 0) then
  !       sp = pspc1
  !     else
  !       sp = pspc2
  !     end if
  !     iy1 = int((brdy2-brdy1)/sp/2)
  !     iz1 = int((brdz2-brdz1)/sp/2)
  !     ! print*,(-by -iy1), (iy1 +by), by, iy1
  !     ! print*, '------'
  !     do j = (-by -iy1), (iy1 +by)
  !       ! print*, j
  !       do k = (-bz -iz1), (iz1 +bz)
  !         isborder = 0
  !         if (i < -ix1) then
  !           isborder = 1
  !           call bx1%append(n)
  !         else if (i < -ix1 + bx) then
  !           call inbx1%append(n)
  !         else if (i <= ix2 - bx) then
  !         else if (i <= ix2) then
  !           call inbx2%append(n)
  !         else
  !           isborder = 1
  !           call bx2%append(n)
  !         end if
  !         if (dim > 1) then
  !           if (j < -iy1) then
  !             isborder = 1
  !             call by1%append(n)
  !           else if (j < -iy1 + by) then
  !             call inby1%append(n)
  !           else if (j <= iy1 - by) then
  !           else if (j <= iy1) then
  !             call inby2%append(n)
  !           else
  !             isborder = 1
  !             call by2%append(n)
  !           end if
  !           if (dim == 3) then
  !             if (k < -iz1) then
  !               isborder = 1
  !               call bz1%append(n)
  !             else if (k < -iz1 + bz) then
  !               call inbz1%append(n)
  !             else if (k <= iz1 - bz) then
  !             else if (k <= iz1) then
  !               call inbz2%append(n)
  !             else
  !               isborder = 1
  !               call bz2%append(n)
  !             end if
  !           end if
  !         end if
  !         ptsz = size(pos, dim=2)
  !         if (ptsz < n) then
  !           call resize(pos, ptsz, ptsz*2)
  !           call resize(ptype, ptsz, ptsz*2)
  !         end if
  !         pos(1,n) = sp*i
  !         pos(2,n) = sp*j
  !         pos(3,n) = sp*k
  !         freenumber = merge(freenumber + 1, freenumber, isborder == 0)
  !         ptype(n) = merge(1, 0, isborder == 0)
  !         n = n + 1
  !       end do
  !     end do
  !     ! read*
  !   end do
  !
  !   n = n - 1
  !
  !   call resize(pos,n,n)
  !   call resize(ptype,n,n)
  !
  !   write(*, "(A, F7.5, A, F7.5)") " # #        actual dx:   x1=", pspc1, "   x2=", pspc2
  !   write(*, "(A, I16, A, I16, A)") " # #         p.number:  ", n, &
  !             " total   ", freenumber, " real"
  !   print *, '# #       border x:', bx1%llen(), bx2%llen()
  !   print *, '# # inner border x:', inbx1%llen(), inbx2%llen()
  !   print *, '# #       border y:', by1%llen(), by2%llen()
  !   print *, '# # inner border y:', inby1%llen(), inby2%llen()
  !   print *, '# #       border z:', bz1%llen(), bz2%llen()
  !   print *, '# # inner border z:', inbz1%llen(), inbz2%llen()
  !
  !   call set_particles_numbers(n, abs(nb))
  !
  !   call setBorder(11, bx1%toarr())
  !   call setBorder(12, bx2%toarr())
  !   call setBorder(21, by1%toarr())
  !   call setBorder(22, by2%toarr())
  !   call setBorder(31, bz1%toarr())
  !   call setBorder(32, bz2%toarr())
  !
  !   call setBorderInside(11, inbx1%toarr())
  !   call setBorderInside(12, inbx2%toarr())
  !   call setBorderInside(21, inby1%toarr())
  !   call setBorderInside(22, inby2%toarr())
  !   call setBorderInside(31, inbz1%toarr())
  !   call setBorderInside(32, inbz2%toarr())
  ! end subroutine
end module
