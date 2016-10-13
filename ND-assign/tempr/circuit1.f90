module tempr_circuit1
  use omp_lib
  use kernel

  implicit none

  public :: tempr_c1

  private
contains

  subroutine tempr_c1(n, pos, mas, sk, sln, den, om)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n,3), mas(n), sk
    real, intent(out)   :: sln(n), den(n), om(n)
    real                :: w, dwdh, r(3), dr, dfdh, fh, hn
    real                :: allowerror, slnint(n), resid(n)
    integer             :: i, j, dim, maxiter

    call get_dim(dim)
    ! print *, 'Dim in circuit1: ', dim
    ! read *
    allowerror = 0.0001
    slnint(:) = sln(:)
    resid(:)  = 1.

    do while (maxval(resid, mask=(resid>0)).gt.allowerror)
      !$OMP PARALLEL
      !$OMP DO PRIVATE(r, dr, dwdh, w, dfdh, fh, hn, maxiter)
      do i = 1, n
        if (resid(i).gt.allowerror) then
          den(i) = 0.
          om(i) = 0.
          j = 0
          if (i.ne.j) then
            do j = 1, n
              r(:) = pos(i,:) - pos(j,:)
              dr = sqrt(dot_product(r(:),r(:)))
              if (dr < 2. * slnint(i)) then
                call get_dw_dh(dr, slnint(i), dwdh)
                call get_w(dr, slnint(i), w)
                den(i) = den(i) + mas(j) * w
                om(i) = om(i) + mas(j) * dwdh
              end if
            end do
          end if
          om(i) = 1. - om(i) * (- slnint(i) / (3 * dim * den(i)))
          dfdh = - 3 * dim * den(i) * om(i) / slnint(i)
          fh  = mas(i) * (sk / slnint(i)) ** dim - den(i)
          hn = slnint(i) - fh / dfdh
          resid(i) = abs(hn - slnint(i)) / sln(i)
          slnint(i) = hn
        end if
      end do
      !$OMP END DO
      !$OMP END PARALLEL
    end do
    sln(:) = slnint(:)
  end subroutine tempr_c1

end module tempr_circuit1
