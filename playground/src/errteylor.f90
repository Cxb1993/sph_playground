module errteylor
  use kernel
  use neighboursearch, only:  getneighbours,&
                              isInitialized,&
                              findneighbours,&
                              findneighboursonce
  implicit none

  public :: laplace, setStepsize

  private
  integer, save :: stepsize = 1

contains
  subroutine setStepsize(i)
    integer, intent(in) :: i
    stepsize = i
  end subroutine setStepsize

  subroutine laplace(idx, pos, mas, den, h, chi)
    real, allocatable, intent(in)    :: mas(:), den(:), pos(:,:), h(:)
    integer, intent(in)              :: idx
    real, intent(inout) :: chi(9)

    integer, allocatable :: nlist(:)
    integer              :: i, j, l, kd, ininb!, n
    real                 :: n2w, r(3), r11, r22, r33, r12, r13, r23, kr

    call get_krad(kr)
    call get_dim(kd)
    call isInitialized(ininb)
    ! n = size(ptype)
    chi(1:9) = 0.
    if (ininb == 0) then
      call findneighboursonce(idx, pos, h, nlist)
    end if
    !<$omp parallel do default(none) &
    !<$omp private(n2w,r,r11,r22,r33,r12,r23,r13,j,i,nlist) &
    !<$omp shared(pos,h,mas,den,kr,kd,n,chi,stepsize,ptype)
    ! do i = 1,n,stepsize
    !   if (ptype(i) /= 0) then
        ! call getneighbours(idx, nlist)
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

          chi(1) = chi(1) + 0.5 * mas(j)/den(j) * r11 * n2w
          if (kd > 1) then
            chi(2) = chi(2) + 0.5 * mas(j)/den(j) * r12 * n2w
            chi(5) = chi(5) + 0.5 * mas(j)/den(j) * r22 * n2w
            chi(4) = chi(4) + 0.5 * mas(j)/den(j) * r12 * n2w
            if (kd == 3) then
              chi(3) = chi(3) + 0.5 * mas(j)/den(j) * r13 * n2w
              chi(6) = chi(6) + 0.5 * mas(j)/den(j) * r23 * n2w
              chi(9) = chi(9) + 0.5 * mas(j)/den(j) * r33 * n2w
              chi(8) = chi(8) + 0.5 * mas(j)/den(j) * r23 * n2w
              chi(7) = chi(7) + 0.5 * mas(j)/den(j) * r13 * n2w
            end if
          end if
          ! chi(1,i) = chi(1,i) + 0.5 * mas(j)/den(j) * r11 * n2w
          ! if (kd > 1) then
          !   chi(2,i) = chi(2,i) + 0.5 * mas(j)/den(j) * r12 * n2w
          !   chi(5,i) = chi(5,i) + 0.5 * mas(j)/den(j) * r22 * n2w
          !   chi(4,i) = chi(4,i) + 0.5 * mas(j)/den(j) * r12 * n2w
          !   if (kd == 3) then
          !     chi(3,i) = chi(3,i) + 0.5 * mas(j)/den(j) * r13 * n2w
          !     chi(6,i) = chi(6,i) + 0.5 * mas(j)/den(j) * r23 * n2w
          !     chi(9,i) = chi(9,i) + 0.5 * mas(j)/den(j) * r33 * n2w
          !     chi(8,i) = chi(8,i) + 0.5 * mas(j)/den(j) * r23 * n2w
          !     chi(7,i) = chi(7,i) + 0.5 * mas(j)/den(j) * r13 * n2w
          !   end if
          ! end if
        end do
        ! print *, chi(1:3)
        ! print *, chi(4:6)
        ! print *, chi(7:9)
        ! read *
    !   end if
    ! end do
    !<$omp end parallel do
    !
    ! ? ! Second derivatives term is for particle 'a', but ex is sum of all particle
    ! ? ! Need to divide on number of particles to get an average term
    !
  end subroutine laplace
end module errteylor
