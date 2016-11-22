module bias
  use kernel

  implicit none

  public :: ex22

  private

contains
  subroutine ex22(n, mas, den, pos, h, ex)
    integer, intent(in) :: n
    real, intent(in)    :: mas(n), den(n), pos(3,n), h(n)
    real, intent(out)   :: ex(n)
    integer             :: i, j
    real                :: n2w, r(3), r2, kr

    call get_krad(kr)

    !$omp parallel do default(none) &
    !$omp private(n2w, r, r2, j, i) &
    !$omp shared(pos, h, mas, den, kr, n, ex)
    do i = 1,n
      ex(i) = 0.
      do j = 1,n
        if (i /= j) then
          r(:) = pos(:,j) - pos(:,i)
          r2 = dot_product(r(:), r(:))
          if (r2 < (kr * h(i))**2) then
            call get_n2w(r, h(i), n2w)
            ex(i) = ex(i) + 0.5 * mas(j)/den(j) * r2 * n2w
          end if
        end if
      end do
    end do
    !$omp end parallel do
    !
    ! ? ! Second derivatives term is for particle 'a', but ex is sum of all particle
    ! ? ! Need to divide on number of particles to get an average term
    !
  end subroutine ex22
end module bias
