module circuit1
  use omp_lib
  use timing,          only: addTime
  use kernel
  use neighboursearch, only: getneighbours,&
                             getNeibListL1,&
                             getNeibListL2

  implicit none

  public :: c1_init, c1, c1a, setStepsize

  private
  save
    integer :: stepsize = 1
    real, allocatable :: slnint(:), resid(:)
    integer(8) :: start=0, finish=0

contains
  subroutine c1_init(n)
    integer, intent(in) :: n
    allocate(slnint(n))
    allocate(resid(n))
  end subroutine c1_init

  subroutine setStepsize(i)
    integer, intent(in) :: i
    stepsize = i
  end subroutine setStepsize

  subroutine c1(ptype, pos, mas, vel, sk, h, den, om, dfdx)
    real, allocatable, intent(in)    :: pos(:,:), mas(:), vel(:,:)
    real, allocatable, intent(inout) :: h(:), den(:), om(:), dfdx(:,:,:)
    integer, allocatable, intent(in) :: ptype(:)
    real, intent(in)     :: sk
    real                 :: w, dwdh, r(3), dr, r2, dfdh, fh, hn, vba(3), nw(3)
    real                 :: allowerror
    integer              :: n, ni, nj, i, j, la, lb, dim, iter, ktp
    integer(8)           :: t0, tneib
    integer, allocatable :: nlista(:), nlistb(:)
    call system_clock(start)

    n = size(ptype)

    call get_dim(dim)
    call get_kerntype(ktp)

    if ( ktp == 3 ) then
      call getNeibListL2(nlista)
    else
      call getNeibListL1(nlista)
    end if

    allowerror = 1e-8
    slnint(:) = h(:)
    resid(:)  = 1.
    iter = 0
    tneib = 0.

    do while ((maxval(resid, mask=(resid>0)) > allowerror) .and. (iter < 100))
      iter = iter + 1
      !$omp parallel do default(none)&
      !$omp private(r, dr, dwdh, w, dfdh, fh, hn, j, i, la, lb, r2, t0, nlistb)&
      !$omp private(ni, nj, nw, vba)&
      !$omp shared(ptype, resid, allowerror, n, pos, mas, dim, sk, h, ktp)&
      !$omp shared(nlista, den, om, slnint, dfdx, vel)&
      !$omp reduction(+:tneib)
      do la = 1, size(nlista)
        i = nlista(la)
        if (resid(i) > allowerror) then
          den(i)  = 0.
          om(i)   = 0.
          dfdx(:,:,i) = 0.
          call getneighbours(i, pos, h, nlistb, t0)
          tneib = tneib + t0
          do lb = 1, size(nlistb)
            j = nlistb(lb)
            r(:) = pos(:,i) - pos(:,j)
            r2 = dot_product(r(:),r(:))
            dr = sqrt(r2)
            call get_dw_dh(dr, slnint(i), dwdh)
            call get_w(dr, slnint(i), w)
            den(i) = den(i) + mas(j) * w
            om(i) = om(i) + mas(j) * dwdh
            vba(:) = vel(:,j) - vel(:,i)
            call get_nw(r, h(i), nw)
            do ni = 1,dim
              do nj = 1,dim
                dfdx(ni,nj,i) = dfdx(ni,nj,i) + mas(j)*vba(ni)*nw(nj)
              end do
            end do
          end do
          ! ---------------------------------------------
          ! There is no particle itself in neighbour list
          ! ---------------------------------------------
          call get_dw_dh(0., slnint(i), dwdh)
          call get_w(0., slnint(i), w)
          den(i) = den(i) + mas(i) * w
          om(i) = om(i) + mas(i) * dwdh
          ! ---------------------------------------------
          om(i) = 1. - om(i) * (- slnint(i) / (dim * den(i)))
          dfdx(:,:,i) = dfdx(:,:,i) / om(i) /den(i)
          dfdh = - dim * den(i) * om(i) / slnint(i)
          fh  = mas(i) * (sk / slnint(i)) ** dim - den(i)
          hn = slnint(i) - fh / dfdh
          resid(i) = abs(hn - slnint(i)) / h(i)
          slnint(i) = hn
        end if
      end do
      !$omp end parallel do
    end do
    ! do la = 1,size(nlista)
    !   i = nlista(la)
    !   print*, i, den(i)
    ! end do
    h(:) = slnint(:)
    call system_clock(finish)
    call addTime(' circuit1', finish - start - tneib)
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
