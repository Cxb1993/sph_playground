module errteylor
  use omp_lib

  use timing,           only: addTime
  use kernel,           only: n2w, &
                              hessian, &
                              get_krad
  use state,            only: getdim
  use neighboursearch,  only: getneighbours,&
                              isInitialized
  use BC,               only: getSqaureBoxSides
  use printer,          only: AppendLine

  implicit none

  public :: laplace, graddiv, setStepsize, setInfluenceCalc

  private
  save
  integer    :: stepsize = 1
  integer(8) :: start=0, finish=0
  character (len=100) :: kinfname = ''

contains
  subroutine setStepsize(i)
    integer, intent(in) :: i
    stepsize = i
  end subroutine

  subroutine setInfluenceCalc(fname)
    character (len=*), intent(in) :: fname
    kinfname = adjustl(fname)
  end subroutine

  subroutine laplace(pos, mas, den, h, chi)
    real, allocatable, intent(in)    :: mas(:), den(:), pos(:,:), h(:)
    real, intent(inout) :: chi(9)

    integer, allocatable :: nlist(:)
    integer              :: i, j, l, kd, nx, ny, nz, idx
    integer(8)           :: tneib, tprint
    real                 :: n2wa, r(3), r11, r22, r33, r12, r13, r23, kr, t(3),&
                            dsum, osum, dchi(9)
    real, allocatable    :: printres(:)

    allocate(printres(6))

    call system_clock(start)

    ! print*, kinfname

    call get_krad(kr)
    call getdim(kd)
    chi(1:9) = 0.
    dchi(:) = 0.
    t(:) = 0.
    call getSqaureBoxSides(nx, ny, nz)
    if (kd == 1) then
      idx = int(nx/2)
    else if (kd == 2) then
      idx = int(ny/2*(nx+1))
    else if (kd == 3) then
      idx = int(nz/2*(ny*(nx+1)+1))
    end if
    call getneighbours(idx, pos, h, nlist, tneib)
    i = idx

    dsum = 0.
    osum = 0.

    do l = 1,size(nlist)
      j = nlist(l)
      r(:) = pos(:,j) - pos(:,i)
      r11 = r(1)*r(1)
      r22 = r(2)*r(2)
      r33 = r(3)*r(3)
      r12 = r(1)*r(2)
      r13 = r(1)*r(3)
      r23 = r(2)*r(3)
      call n2w(r, h(i), n2wa)

      ! call GradDivW(r, h(i), n2wa)
      ! call get_Hesobian(r, h(i), Hes)
      t(:) = t(:) + mas(j)/den(j) * r(:) * n2wa

      dchi(1) = 0.5 * mas(j)/den(j) * r11 * n2wa ! Hes(1,1)
      if ( kd > 1 ) then
        dchi(2) = 0.5 * mas(j)/den(j) * r12 * n2wa ! Hes(1,2)
        dchi(5) = 0.5 * mas(j)/den(j) * r22 * n2wa ! Hes(2,2)
        dchi(4) = 0.5 * mas(j)/den(j) * r12 * n2wa ! Hes(1,2)
      end if
      if ( kd == 3 ) then
        dchi(3) = 0.5 * mas(j)/den(j) * r13 * n2wa ! Hes(1,3)
        dchi(6) = 0.5 * mas(j)/den(j) * r23 * n2wa ! Hes(2,3)
        dchi(9) = 0.5 * mas(j)/den(j) * r33 * n2wa ! Hes(3,3)
        dchi(8) = 0.5 * mas(j)/den(j) * r23 * n2wa ! Hes(2,3)
        dchi(7) = 0.5 * mas(j)/den(j) * r13 * n2wa ! Hes(1,3)
      end if
      chi(:) = chi(:) + dchi(:)
      if ( kinfname /= '' ) then
        dsum = dchi(1) + dchi(5) + dchi(9)
        osum = dchi(2) + dchi(4) + dchi(3) + dchi(6) + dchi(8) + dchi(7)
        printres(1:3) = r(:)
        printres(4) = sqrt(dot_product(r,r))
        printres(5) = dsum
        printres(6) = osum
        call AppendLine(printres, kinfname, tprint)
        tneib = tneib + tprint
        ! print*, r(:), sqrt(dot_product(r,r)), dsum, osum
      end if

      ! print*, n2w

    end do
    ! print*,
    ! stop
    ! print*, ' t: ', t
    call system_clock(finish)
    call addTime(' teylor err', finish - start - tneib)
  end subroutine laplace

  subroutine graddiv(pos, mas, den, h, chi)
    real, allocatable, intent(in)    :: mas(:), den(:), pos(:,:), h(:)
    real, intent(inout)              :: chi(81)

    integer, allocatable :: nlist(:)
    integer              :: i, j, l, kd, nx, ny, nz, idx, a, b, g, d, ci
    integer(8)           :: tneib, tprint
    real                 :: r(3), kr, t(3), dr, Hes(3,3), m, dchi(81), dsum, osum
    real, allocatable    :: printres(:)

    allocate(printres(6))


    call system_clock(start)

    call get_krad(kr)
    call getdim(kd)
    chi(1:81) = 0.
    dchi(:) = 0.
    t(:) = 0.

    call getSqaureBoxSides(nx, ny, nz)
    if (kd == 1) then
      idx = int(nx/2)
    else if (kd == 2) then
      idx = int(ny/2*(nx+1))
    else if (kd == 3) then
      idx = int(nz/2*(ny*(nx+1)+1))
    end if
    call getneighbours(idx, pos, h, nlist, tneib)
    i = idx
    do l = 1,size(nlist)
      j = nlist(l)
      r(:) = pos(:,i) - pos(:,j)
      call hessian(r, pos(:,i), pos(:,j), h(i), Hes)
      ci = 1
      do a = 1,3
        do b = 1,3
          dr = r(a)*r(b)
          if (a == b) then
            m = .5
          else
            m = 1.
          end if
          do g = 1,3
            do d = 1,3
              dchi(ci) = m * mas(j) / den(j) * dr * Hes(g,d)
              ci = ci + 1
            end do
          end do
        end do
      end do
      chi(:) = chi(:) + dchi(:)
      if ( kinfname /= '' ) then
        dsum = 0.
        osum = 0.
        ci = 1
        do a = 1,3
          do b = 1,3
            do g = 1,3
              do d = 1,3
                if ((( a == g ).and.( b == d )).or.(( a == d ).and.( b == g ))) then
                  dsum = dsum + dchi(ci)
                else
                  osum = osum + dchi(ci)
                end if
                ci = ci + 1
              end do
            end do
          end do
        end do
        printres(1:3) = r(:)
        printres(4) = sqrt(dot_product(r,r))
        printres(5) = dsum
        printres(6) = osum
        call AppendLine(printres, kinfname, tprint)
        tneib = tneib + tprint
      end if
    end do
    call system_clock(finish)
    call addTime(' teylor err', finish - start - tneib)
  end subroutine graddiv
end module
