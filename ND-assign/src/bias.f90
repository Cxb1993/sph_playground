module bias
  use kernel

  implicit none

  public :: ex22

  private

contains
  subroutine ex22(n, mas, den, pos, h, ex)
    integer, intent(in) :: n
    real, intent(in)    :: mas(n), den(n), pos(3,n), h(n)
    real, intent(out)   :: ex
    integer             :: i, j
    character(len=40)   :: kt
    real                :: n2w, r(3), r2, kr

    call get_kerntype(kt)
    call get_krad(kr)

    ex = 0

    !$omp parallel do default(none) &
    !$omp private(n2w, r, r2, j, i) &
    !$omp shared(pos, h, mas, den, kt, kr, n) &
    !$omp reduction(+:ex)
    do i = 1,n
      do j = 1,n
        if (i /= j) then
          r(:) = pos(:,j) - pos(:,i)
          r2 = dot_product(r(:), r(:))
          if (r2 < kr * kr * h(i) * h(i)) then
            if (kt == 'n2w') then
              call get_n2w(sqrt(r2), h(i), n2w)
            else if (kt == 'fab') then
              call get_Fab(r, h(i), n2w)
            else
              print *, 'kernel type not chosen, arg #5'
              stop
            end if
            ex = ex + mas(j)/den(j) * r2 * n2w
          end if
        end if
      end do
    end do
    !$omp end parallel do
    ex = 1./2. * ex / n
  end subroutine ex22
end module bias
