module circuit1
  use omp_lib
  use kernel

  implicit none

  public :: make_c1

  private
contains

  subroutine make_c1(n, pos, mas, sk, sln, den, om)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n,3), mas(n), sk
    real, intent(out)   :: sln(n), den(n), om(n)
    real                :: w, dwdh, r(3), dr, dfdh, fh, hn
    real                :: allowerror, slnint(n), resid(n)
    integer             :: i, j, dim, maxiter

    call get_dim(dim)
    allowerror = 0.0001
    slnint(:) = sln(:)
    resid(:)  = 1.

    do while (maxval(resid, mask=(resid>0)).gt.allowerror)
      ! print *, maxval(resid, mask=(resid>0))
      ! read *
      !$omp parallel do default(none)&
      !$omp private(r, dr, dwdh, w, dfdh, fh, hn, maxiter)&
      !$omp shared(resid, allowerror, den, om, n, pos, slnint, mas, dim, sk, sln)
      do i = 1, n
        if (resid(i) > allowerror) then
          den(i) = 0.
          om(i) = 0.
          do j = 1, n
            if (i.ne.j) then
              r(:) = pos(i,:) - pos(j,:)
              dr = sqrt(dot_product(r(:),r(:)))
              if (dr < 3. * slnint(i)) then
                call get_dw_dh(dr, slnint(i), dwdh)
                ! print *, slnint(i), sln(i)
                ! read *
                call get_w(dr, slnint(i), w)
                den(i) = den(i) + mas(j) * w
                om(i) = om(i) + mas(j) * dwdh
              end if
            end if
          end do
          om(i) = 1. - om(i) * (- slnint(i) / (dim * den(i)))
          dfdh = - dim * den(i) * om(i) / slnint(i)
          fh  = mas(i) * (sk / slnint(i)) ** dim - den(i)
          hn = slnint(i) - fh / dfdh
          resid(i) = abs(hn - slnint(i)) / sln(i)
          slnint(i) = hn
        end if
        ! print *, den(i)
        ! read *
      end do
      !$omp end parallel do
      ! print *, 1, den(1), int(size(den)/2), den(int(size(den)/2)), int(size(den)), den(int(size(den)))
    end do
    sln(:) = slnint(:)
    ! print *, 1, den(1), int(size(den)/2), den(int(size(den)/2)), int(size(den)), den(int(size(den)))
    ! print *, 1, sln(1), int(size(sln)/2), sln(int(size(sln)/2)), int(size(sln)), sln(int(size(sln)))
    ! read *
  end subroutine make_c1

end module circuit1
