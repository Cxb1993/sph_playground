module cltest

contains
  subroutine runcltest()
    use clfortran
    use ISO_C_BINDING
    implicit none

    integer(c_int) :: err
    integer(c_int) :: num_platforms
    integer(c_int) :: i
    integer(c_int64_t) :: device_type
    integer(c_int) :: num_devices
    type(c_ptr), allocatable, target :: platform_ids(:)

    ! Get platform count to allocate an array.
    print *, clGetPlatformIDs(0, C_NULL_PTR, num_platforms)
    print *, 'Num Platforms: ', num_platforms

    ! Allocate an array to hold platforms.
    allocate(platform_ids(num_platforms))

    ! Get platforms IDs.
    err = clGetPlatformIDs(num_platforms, C_LOC(platform_ids), num_platforms)

    ! Loop over all platforms and query devices.
    do i = 1, num_platforms
        ! Iterate over platforms and get number of devices.
        print *, 'Platform: ', i

        ! Get device count.
        device_type = -1        ! ALL, a constant will be added later
        ! print *, clGetDeviceIDs(platform_ids(i), device_type, 0, C_NULL_PTR, num_devices)
        print *, 'Num Devices: ', num_devices
    end do
  end subroutine
end module
