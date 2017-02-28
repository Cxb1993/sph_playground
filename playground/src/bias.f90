module bias
  use kernel
  use neighboursearch, only:  getneighbours,&
                              isInitialized,&
                              findneighbours
  implicit none

  public :: calcDaigonal2ndErrTerms, setStepsize

  private
  integer, save :: stepsize = 1

contains
  subroutine setStepsize(i)
    integer, intent(in) :: i
    stepsize = i
  end subroutine setStepsize

  subroutine calcDaigonal2ndErrTerms(ptype, pos, mas, den, h, chi)
    real, allocatable, intent(in)  :: mas(:), den(:), pos(:,:), h(:)
    integer, allocatable, intent(in) :: ptype(:)
    real, allocatable, intent(inout) :: chi(:,:)

    integer, allocatable :: nlist(:)
    integer              :: i, j, l, kd, n, ininb
    real                 :: n2w(3), r(3), r2, kr

    call get_krad(kr)
    call get_dim(kd)
    call isInitialized(ininb)
    n = size(ptype)
    chi(1:3,1:n) = 0.
    if (ininb == 0) then
      call findneighbours(ptype, pos, h)
    end if
    !$omp parallel do default(none) &
    !$omp private(n2w,r,r2,j,i,nlist) &
    !$omp shared(pos,h,mas,den,kr,kd,n,chi,stepsize,ptype)
    do i = 1,n,stepsize
      if (ptype(i) /= 0) then
        call getneighbours(i,nlist)
        do l = 1,size(nlist)
          j = nlist(l)
          r(:) = pos(:,j) - pos(:,i)
          r2 = dot_product(r(:), r(:))
          call GradDivW(r, h(i), n2w)
          ! call get_n2iw(r, h(i), n2w(1), 1)
          ! if (kd > 1) then
          !   call get_n2iw(r, h(i), n2w(2), 2)
          !   if (kd == 3) then
          !     call get_n2iw(r, h(i), n2w(3), 3)
          !   end if
          ! end if
          chi(:,i) = chi(:,i) + 0.5 * mas(j)/den(j) * r2 * n2w(:)
        end do
      end if
    end do
    !$omp end parallel do
    !
    ! ? ! Second derivatives term is for particle 'a', but ex is sum of all particle
    ! ? ! Need to divide on number of particles to get an average term
    !
  end subroutine calcDaigonal2ndErrTerms
end module bias
