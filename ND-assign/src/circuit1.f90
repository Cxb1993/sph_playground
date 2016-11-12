module circuit1
  use omp_lib
  use kernel

  implicit none

  public :: c1, c1a

  private
contains

  subroutine c1(n, pos, mas, sk, sln, den, om)
    integer, intent(in) :: n
    real, intent(in)    :: pos(3,n), mas(n), sk
    real, intent(out)   :: sln(n), den(n), om(n)
    real                :: w, dwdh, r(3), dr, r2, dfdh, fh, hn, kr
    real                :: allowerror, slnint(n), resid(n)
    integer             :: i, j, dim, maxiter

    call get_dim(dim)
    call get_krad(kr)

    allowerror = 1e-8
    slnint(:) = sln(:)
    resid(:)  = 1.
    do while (maxval(resid, mask=(resid>0)).gt.allowerror)
      !$omp parallel do default(none)&
      !$omp private(r, dr, dwdh, w, dfdh, fh, hn, maxiter, j, i, r2)&
      !$omp shared(resid, allowerror, den, om, n, pos, slnint, mas, dim, sk, sln, kr)
      do i = 1, n
        if (resid(i) > allowerror) then
          den(i) = 0.
          om(i) = 0.
          do j = 1, n
            r(:) = pos(:,i) - pos(:,j)
            r2 = dot_product(r(:),r(:))
            if (r2 <= (kr * slnint(i))**2) then
              dr = sqrt(r2)
              call get_dw_dh(dr, slnint(i), dwdh)
              call get_w(dr, slnint(i), w)
              den(i) = den(i) + mas(j) * w
              ! print *, 'in sum: ', mas(j) * w, 'mas: ', mas(j), 'kernel: ', w
              om(i) = om(i) + mas(j) * dwdh
            end if
          end do
          om(i) = 1. - om(i) * (- slnint(i) / (dim * den(i)))
          dfdh = - dim * den(i) * om(i) / slnint(i)
          fh  = mas(i) * (sk / slnint(i)) ** dim - den(i)
          hn = slnint(i) - fh / dfdh
          resid(i) = abs(hn - slnint(i)) / sln(i)
          slnint(i) = hn
        end if
      end do
      !$omp end parallel do
    end do
    sln(:) = slnint(:)
  end subroutine c1

! Direct density summation
  subroutine c1a(n, pos, mas, sk, sln, den)
    integer, intent(in) :: n
    real, intent(in)    :: pos(3,n), mas(n), sk
    real, intent(out)   :: den(n), sln(n)
    real                :: w, r(3), dr
    integer             :: i, j

    do i = 1, n
      den(i) = 0.
      do j = 1, n
        r(:) = pos(:,i) - pos(:,j)
        dr = sqrt(dot_product(r(:),r(:)))
        if (dr <= 2. * sln(i)) then
          call get_w(dr, sln(i), w)
          den(i) = den(i) + mas(j) * w
        endif
      end do
      sln(i) = sk * (mas(i) / den(i))
    end do
  end subroutine c1a
end module circuit1
