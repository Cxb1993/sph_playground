module errteylor
  use timing,          only: addTime
  use kernel
  use neighboursearch, only:  getneighbours,&
                              isInitialized,&
                              findneighbours,&
                              getneighbours
  use BC,              only: getSqaureBoxSides
  implicit none

  public :: laplace, graddiv, setStepsize

  private
  integer, save :: stepsize = 1
  integer(8) :: start=0, finish=0

contains
  subroutine setStepsize(i)
    integer, intent(in) :: i
    stepsize = i
  end subroutine setStepsize

  subroutine laplace(pos, mas, den, h, chi)
    real, allocatable, intent(in)    :: mas(:), den(:), pos(:,:), h(:)
    real, intent(inout) :: chi(9)

    integer, allocatable :: nlist(:)
    integer              :: i, j, l, kd, nx, ny, nz, idx
    integer(8)           :: tneib
    real                 :: n2w, r(3), r11, r22, r33, r12, r13, r23, kr, t(3)!, n2wa(3)!, Hes(3,3)

    call system_clock(start)

    call get_krad(kr)
    call get_dim(kd)
    chi(1:9) = 0.
    t(:) = 0.

    call getSqaureBoxSides(nx, ny, nz)
    if (kd == 1) then
      idx = int(nx/2)
    else if (kd == 2) then
      idx = int(ny/2*(nx+1))
    else if (kd == 3) then
      idx = int(nz/2*(ny*(nx+1)+1))
    end if
    call getneighbours(idx, pos, h, nlist, tneib)
    i = idx
    do l = 1,size(nlist)
      j = nlist(l)
      r(:) = pos(:,j) - pos(:,i)
      r11 = r(1)*r(1)
      r22 = r(2)*r(2)
      r33 = r(3)*r(3)
      r12 = r(1)*r(2)
      r13 = r(1)*r(3)
      r23 = r(2)*r(3)
      call get_n2w(r, h(i), n2w)
      ! call GradDivW(r, h(i), n2wa)
      ! call get_Hesobian(r, h(i), Hes)
      t(:) = t(:) + mas(j)/den(j) * r(:) * n2w

      chi(1) = chi(1) + 0.5 * mas(j)/den(j) * r11 * n2w ! Hes(1,1)
      if (kd > 1) then
        chi(2) = chi(2) + 0.5 * mas(j)/den(j) * r12 * n2w ! Hes(1,2)
        chi(5) = chi(5) + 0.5 * mas(j)/den(j) * r22 * n2w ! Hes(2,2)
        chi(4) = chi(4) + 0.5 * mas(j)/den(j) * r12 * n2w ! Hes(1,2)
        if (kd == 3) then
          chi(3) = chi(3) + 0.5 * mas(j)/den(j) * r13 * n2w ! Hes(1,3)
          chi(6) = chi(6) + 0.5 * mas(j)/den(j) * r23 * n2w ! Hes(2,3)
          chi(9) = chi(9) + 0.5 * mas(j)/den(j) * r33 * n2w ! Hes(3,3)
          chi(8) = chi(8) + 0.5 * mas(j)/den(j) * r23 * n2w ! Hes(2,3)
          chi(7) = chi(7) + 0.5 * mas(j)/den(j) * r13 * n2w ! Hes(1,3)
        end if
      end if
    end do
    ! print*, ' t: ', t
    call system_clock(finish)
    call addTime(' teylor err', finish - start - tneib)
  end subroutine laplace

  subroutine graddiv(pos, mas, den, h, chi)
    real, allocatable, intent(in)    :: mas(:), den(:), pos(:,:), h(:)
    real, intent(inout)              :: chi(81)

    integer, allocatable :: nlist(:)
    integer              :: i, j, l, kd, nx, ny, nz, idx, a, b, g, d, ci
    integer(8)           :: tneib
    real                 :: r(3), kr, t(3), dr, Hes(3,3), m

    call system_clock(start)

    call get_krad(kr)
    call get_dim(kd)
    chi(1:81) = 0.
    t(:) = 0.

    call getSqaureBoxSides(nx, ny, nz)
    if (kd == 1) then
      idx = int(nx/2)
    else if (kd == 2) then
      idx = int(ny/2*(nx+1))
    else if (kd == 3) then
      idx = int(nz/2*(ny*(nx+1)+1))
    end if
    call getneighbours(idx, pos, h, nlist, tneib)
    i = idx
    do l = 1,size(nlist)
      j = nlist(l)
      r(:) = pos(:,j) - pos(:,i)
      call get_hessian(r, h(i), Hes)
      ci = 1
      do a = 1,3
        do b = 1,3
          dr = r(a)*r(b)
          if (a == b) then
            m = .5
          else
            m = 1.
          end if
          do g = 1,3
            do d = 1,3
              chi(ci) = chi(ci) + m * mas(j) / den(j) * dr * Hes(g,d)
              ci = ci + 1
            end do
          end do
        end do
      end do
    end do
    ! ci = 1
    ! do a = 1,3
    !   do b = 1,3
    !     do g = 1,3
    !       do d = 1,3
    !         write(*, "(I2, A, I1, A, I1, A, I1, I1, A, F10.7)") ci, " # 1/2(dr", a,"*dr", b, ")*W", g, d, "=", chi(ci)
    !         ci = ci + 1
    !       end do
    !     end do
    !   end do
    ! end do
    call system_clock(finish)
    call addTime(' teylor err', finish - start - tneib)
  end subroutine graddiv
end module errteylor
