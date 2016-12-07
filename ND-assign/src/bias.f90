module bias
  use kernel

  implicit none

  public :: get_chi

  private

contains
  subroutine get_chi(n, mas, den, pos, h, chi)
    integer, intent(in) :: n
    real, intent(in)    :: mas(n), den(n), pos(3,n), h(n)
    real, intent(out)   :: chi(3,n)
    integer             :: i, j, kd
    real                :: n2w(3), r(3), r2, kr

    call get_krad(kr)
    call get_dim(kd)

    !$omp parallel do default(none) &
    !$omp private(n2w, r, r2, j, i) &
    !$omp shared(pos, h, mas, den, kr, kd, n, chi)
    do i = 1,n
      chi(:,i) = 0.
      do j = 1,n
        if (i /= j) then
          r(:) = pos(:,j) - pos(:,i)
          r2 = dot_product(r(:), r(:))
          if (r2 < (kr * h(i))**2) then
            call get_n2iw(r, h(i), n2w(1), 1)
            if (kd > 1) then
              call get_n2iw(r, h(i), n2w(2), 2)
              if (kd == 3) then
                call get_n2iw(r, h(i), n2w(3), 3)
              end if
            end if
            chi(:,i) = chi(:,i) + 0.5 * mas(j)/den(j) * r2 * n2w(:)
          end if
        end if
      end do
    end do
    !$omp end parallel do
    !
    ! ? ! Second derivatives term is for particle 'a', but ex is sum of all particle
    ! ? ! Need to divide on number of particles to get an average term
    !
  end subroutine get_chi
end module bias
