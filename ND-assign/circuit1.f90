module circuit1_mod
  use omp_lib
  use kernel

  implicit none

  public :: circuit1

  private
contains

  subroutine circuit1(n, pos, mas, sk, sln, den, om)
    integer, intent(in) :: n
    real, intent(in)    :: pos(n,3), mas(n), sk
    real, intent(out)   :: sln(n), den(n), om(n)
    real                :: w, dwdh, r(3), dr, dnrfh, nrfh, h0, hp, hn
    integer             :: i, j, dim, itern

    call get_dim(dim)

    !$OMP PARALLEL
    !$OMP DO PRIVATE(r, dr, dwdh, w, dnrfh, nrfh, h0, hp, hn, itern)
    do i = 1, n
      den(i) = 0.
      om(i) = 0.
      j = 0
      if (i.ne.j) then
        do j = 1, n
          r(:) = pos(i,:) - pos(j,:)
          dr = sqrt(dot_product(r(:),r(:)))
          if (dr < 2. * sln(i)) then
            call get_dw_dh(dr, sln(i), dwdh)
            call get_w(dr, sln(i), w)
            den(i) = den(i) + mas(j) * w
            om(i) = om(i) + mas(j) * dwdh
          end if
        end do
      end if
      om(i) = 1. - om(i) * (- sln(i) / (3 * dim * den(i)))
      h0 = sln(i)
      hp = 1000 * h0
      hn = h0
      itern = 50
      do while ((abs(hn - hp) / h0 > 0.0001).and.(itern.gt.0))
        hp = hn
        dnrfh = - 3 * dim * den(i) * om(i) / hp
        nrfh  = mas(i) * (sk / hp) ** dim - den(i)
        hn = hp - nrfh / dnrfh
        itern = itern - 1
        if (itern.eq.0) then
          print *, 'More then 50 iterations in NR scheme. i:', i, 'h_new', hn
          print *, pos(i,:)
          print *, hp, nrfh, dnrfh
          stop
        end if
      end do
      sln(i) = hn
    end do
    !$OMP END DO
    !$OMP END PARALLEL
  end subroutine circuit1

end module circuit1_mod
