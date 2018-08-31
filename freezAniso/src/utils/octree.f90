module octree

  implicit none

  public :: create

  private

  save
    integer, parameter :: keylen = 16

    type charintmap
      character(len=keylen), allocatable :: key(:)
      integer(8), allocatable  :: val(:)
    end type

    type charintlib
      character(len=keylen), allocatable :: name(:)
      type(charintmap), allocatable :: data(:)
    end type

    type charcharmap
      character(len=keylen), allocatable :: key(:)
      character(len=keylen), allocatable :: val(:)
    end type

    type charcharlib
      character(len=keylen), allocatable :: name(:)
      type(charcharmap), allocatable :: data(:)
    end type

    type(charcharlib) :: ccl
    type(charintlib)  :: cil

  contains

    subroutine expandcclib()
      type(charcharmap), allocatable :: tmpccl(:)
      integer :: n, i

      if (allocated(ccl%data)) then
        n = size(ccl%data)
        allocate(tmpccl(n))
        do i = 1,n
          tmpccl(i) = ccl%data(i)
        end do
        deallocate(ccl%data)
        allocate(ccl%data(n+1))
        do i = 1,n
          ccl%data(i) = tmpccl(i)
        end do
        deallocate(tmpccl)
      else
        allocate(ccl%data(1))
      end if
    end subroutine

    subroutine expandcilib()
      type(charintmap), allocatable :: tmpcil(:)
      integer :: n, i

      if (allocated(cil%data)) then
        n = size(cil%data)
        allocate(tmpcil(n))
        do i = 1,n
          tmpcil(i) = cil%data(i)
        end do
        deallocate(cil%data)
        allocate(cil%data(n+1))
        do i = 1,n
          cil%data(i) = tmpcil(i)
        end do
        deallocate(tmpcil)
      else
        allocate(cil%data(1))
      end if
    end subroutine

    subroutine makemap(inname, inmaptype)
      character(len=*), intent(in) :: inname, inmaptype
      integer :: n, i

      select case(inmaptype)
      case('charchar')
        if ( allocated(ccl%name) ) then
          n = size(ccl%name)
          do i = 1,n
            if ( ccl%name(i) == inname ) then
              print*, 'Such map already has been allocated: ', ccl%name(i)
              stop
            end if
          end do
          call resize(ccl%name, keylen, n, n+1)
          call expandcclib()
          ccl%name(n+1) = inname
        else
          allocate(ccl%name(1))
          call expandcclib()
          ccl%name(1) = inname
        end if
      case('charint')
        if ( allocated(cil%name) ) then
          n = size(cil%name)
          do i = 1,n
            if ( cil%name(i) == inname ) then
              print*, 'Such map already has been allocated: ', cil%name(i)
              stop
            end if
          end do
          call resize(cil%name, keylen, n, n+1)
          call expandcclib()
          cil%name(n+1) = inname
        else
          allocate(cil%name(1))
          call expandcilib()
          cil%name(1) = inname
        end if
      case default
        print *, 'There is no such mapionary type: ', inmaptype
        stop
      end select
    end subroutine

    subroutine appendmapcc(mapname, inkey, inval)
      character(len=*), intent(in) :: mapname
      character(len=*), intent(in) :: inkey
      character(len=*), intent(in) :: inval

      integer :: nl, il, nm, im, fk

      fk = 0

      if (.not. allocated(ccl%name) ) then
        call makemap(mapname, 'charchar')
      end if

      nl = size(ccl%name)

      do il = 1,nl
        if ( ccl%name(il) == mapname ) then
          if (allocated(ccl%data(il)%key)) then
            nm = size(ccl%data(il)%key)
            do im = 1,nm
              if ( ccl%data(il)%key(im) == inkey ) then
                fk = 1
                ccl%data(il)%val(im) = inval
              end if
            end do
            if ( fk == 0 ) then
              call resize(ccl%data(il)%key, keylen, nm, nm + 1)
              call resize(ccl%data(il)%val, keylen, nm, nm + 1)
              ccl%data(il)%key(nm + 1) = inkey
              ccl%data(il)%val(nm + 1) = inval
            end if
          else
            allocate(ccl%data(il)%key(1))
            ccl%data(il)%key(1) = inkey
            allocate(ccl%data(il)%val(1))
            ccl%data(il)%val(1) = inval
          end if
        end if
      end do
    end subroutine

    subroutine appendmapci(mapname, inkey, inval)
      character(len=*), intent(in)  :: mapname
      character(len=*), intent(in)  :: inkey
      integer, intent(in)           :: inval

      integer :: nl, il, nm, im, fk

      fk = 0

      if (.not. allocated(cil%name) ) then
        call makemap(mapname, 'charint')
      end if

      nl = size(cil%name)

      do il = 1,nl
        if ( cil%name(il) == mapname ) then
          if (allocated(cil%data(il)%key)) then
            nm = size(cil%data(il)%key)
            do im = 1,nm
              if ( cil%data(il)%key(im) == inkey ) then
                fk = 1
                cil%data(il)%val(im) = inval
              end if
            end do
            if ( fk == 0 ) then
              call resize(cil%data(il)%key, keylen, nm, nm + 1)
              call resize(cil%data(il)%val, nm, nm + 1)
              cil%data(il)%key(nm + 1) = inkey
              cil%data(il)%val(nm + 1) = inval
            end if
          else
            allocate(cil%data(il)%key(1))
            cil%data(il)%key(1) = inkey
            allocate(cil%data(il)%val(1))
            cil%data(il)%val(1) = inval
          end if
        end if
      end do
    end subroutine

    subroutine destroy()
      integer :: n, i

      if (allocated(ccl%name)) then
        n = size(ccl%name)
        do i = 1,n
          ! print*, ccl%name(i)
          if (allocated(ccl%data(i)%key)) then
            ! print*, ccl%data(i)%key
            deallocate(ccl%data(i)%key)
          end if
          if (allocated(ccl%data(i)%val)) then
            ! print*, ccl%data(i)%val
            deallocate(ccl%data(i)%val)
          end if
        end do
        deallocate(ccl%data)
        deallocate(ccl%name)
      end if

      if (allocated(cil%name)) then
        n = size(cil%name)
        do i = 1,n
          if (allocated(cil%data(i)%key)) then
            deallocate(cil%data(i)%key)
          end if
          if (allocated(cil%data(i)%val)) then
            deallocate(cil%data(i)%val)
          end if
        end do
        deallocate(cil%data)
        deallocate(cil%name)
      end if
    end subroutine
end module
