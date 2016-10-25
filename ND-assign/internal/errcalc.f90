module err_calc
  use omp_lib
  use BC

  implicit none

  public :: err_T0sxsyet

  private
  real, parameter     :: pi = 4.*atan(1.)

contains
  subroutine err_T0sxsyet(n, pos, num, t, err)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n,3), num(n), t
    real, intent(out)   :: err(n)
    integer             :: i

    !$OMP PARALLEL
    !$OMP DO
    do i=1,n
      ! err(i) = abs(sin(pi*(pos(i,1)+1.)/2.) * sin(pi*(pos(i,2)+1.)/2.) * exp(-2.*(pi/2.)**2) * 0.1 * t - num(i))
      err(i) = abs(sin(pi*(pos(i,1)+1.)/2.) * exp(-2.*(pi/2.)**2) * 0.1 * t - num(i))
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine err_T0sxsyet
end module err_calc
