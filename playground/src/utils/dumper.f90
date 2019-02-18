module dumper
  use const
  use timing, only: addTime
  use state,  only: setState

  implicit none

  public :: dump, restore, clean

  private
    integer(8) :: start=0, finish=0
  contains

    subroutine dump(store, t)
      use state, only: setStateVal, getStateVal, getState

      real, allocatable, intent(in) :: store(:,:)
      real, intent(in) :: t

      real :: ifile
      real, allocatable :: state(:)
      character (len=40)  :: cfile
      allocate(state(ec_total))
      call system_clock(start)

      call setStateVal(ec_time, t)
      call getStateVal(ec_lastprint, ifile)
      call getState(state)

      call CreateNewDump("output/new_dump_full.h5")
      call PutData1D("output/new_dump_full.h5",&
        "state", ec_total, state)
      call PutData2D("output/new_dump_full.h5",&
        "store", es_total, int(state(ec_realpn))+int(state(ec_fixedpn)), store)

      call rename("output/new_dump_full.h5", "output/dump_full.h5")
      write(cfile, "(a,i5.5)") "cp output/dump_full.h5 output/dump_", int(ifile)
      call system(cfile)

      call system_clock(finish)
      call addTime('dumper', finish - start)
    end subroutine dump

    subroutine restore(store)
      real, allocatable, intent(inout) :: store(:,:)

      real, allocatable :: state(:)

      call system_clock(start)

      print*, "#  # ---   ---   ---   ---   ---   ---   ---   ---   ---"

      call GetData1D("output/dump_full.h5", "state", state)

      call setState(state)

      write(*,blockFormatFltSci) " #  #", " restored from t = ", state(ec_time)
      write(*,blockFormatFltSci) " #  #", "             dump #", state(ec_lastprint)

      call GetData2D("output/dump_full.h5", "store", store)

      print*, "#  # ---   ---   ---   ---   ---   ---   ---   ---   ---"

      call system_clock(finish)
      call addTime('dumper', finish - start)
    end subroutine restore

    subroutine clean()
      integer :: nu, snu
      open(newunit=nu, iostat=snu, file="output/new_dump_full.h5", status='old')
      if (snu == 0) close(nu, status='delete')
    end subroutine clean

    subroutine CreateNewDump(filename)
      use HDF5
      character(len=*), intent(in) :: filename
      integer(HID_T) :: file_id
      integer :: error

      call h5open_f(error)
      call h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)
    end subroutine

    subroutine PutData1D(file, dspace, dim, data)
      use HDF5
      character(len=*), intent(in) :: file, dspace
      real, allocatable, intent(in) :: data(:)
      INTEGER(HSIZE_T) :: ddim(1), maxdims(1)
      integer, intent(in) :: dim

      integer(HID_T) :: file_id, dspace_id, dset_id, crp_id
      integer :: error

      ! open handler
      call h5open_f(error)
      ! open file
      call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error)
      !prepare dims
      ddim(1) = dim
      maxdims(1) = H5S_UNLIMITED_F
      ! create extendible dataspace
      call h5screate_simple_f(1, ddim, dspace_id, error, maxdims)
      ! get dataset creation property identifier
      call h5pcreate_f(H5P_DATASET_CREATE_F, crp_id, error)
      ! enable data to be in chank
      call h5pset_chunk_f(crp_id, 1, ddim, error)
      ! create double dataspace
      call h5dcreate_f(file_id, dspace, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error, crp_id)
      ! write data
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, ddim, error)
      ! close all the handlers
      call h5pclose_f(crp_id, error)
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)
    end subroutine
    subroutine PutData2D(file, dspace, dim1, dim2, data)
      use HDF5
      character(len=*), intent(in) :: file, dspace
      real, allocatable, intent(in) :: data(:,:)
      INTEGER(HSIZE_T) :: ddim(2), maxdims(2)
      integer, intent(in) :: dim1, dim2

      integer(HID_T) :: file_id, dspace_id, dset_id, crp_id
      integer :: error

      call h5open_f(error)
      call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error)
      ddim = [dim1, dim2]
      maxdims = [H5S_UNLIMITED_F, H5S_UNLIMITED_F]
      call h5screate_simple_f(2, ddim, dspace_id, error, maxdims)
      ! get dataset creation property identifier
      call h5pcreate_f(H5P_DATASET_CREATE_F, crp_id, error)
      ! enable data to be in chank
      call h5pset_chunk_f(crp_id, 2, ddim, error)
      ! create double dataspace
      call h5dcreate_f(file_id, dspace, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error, crp_id)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, ddim, error)

      call h5pclose_f(crp_id, error)
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)
    end subroutine

    subroutine GetData1D(file, dsetname, databuf)
      use HDF5
      character(len=*), intent(in) :: file, dsetname
      real, allocatable, intent(inout) :: databuf(:)

      integer(HSIZE_T) :: dimsr(1), maxdimsr(1)
      integer(HID_T) :: file_id, dspace_id, dset_id, crp_id, mem_id
      integer :: error, rankr

      !Init HDF
      call h5open_f(error)
      !Open the file.
      call h5fopen_f (file, H5F_ACC_RDONLY_F, file_id, error)
      !Open the  dataset.
      call h5dopen_f(file_id, dsetname, dset_id, error)
      !Get dataset's dataspace handle.
      call h5dget_space_f(dset_id, dspace_id, error)
      !Get dataspace's rank.
      call h5sget_simple_extent_ndims_f(dspace_id, rankr, error)
      !Get dataspace's dimensions.
      call h5sget_simple_extent_dims_f(dspace_id, dimsr, maxdimsr, error)
      !Get creation property list.
      CALL h5dget_create_plist_f(dset_id, crp_id, error)
      !create memory dataspace
      call h5screate_simple_f(rankr, dimsr, mem_id, error)
      !Read data
      allocate(databuf(dimsr(1)))
      call H5dread_f(dset_id, H5T_NATIVE_DOUBLE, databuf, dimsr, &
           error, mem_id, dspace_id)
      !Close the objects that were opened.
      call h5sclose_f(dspace_id, error)
      call h5sclose_f(mem_id, error)
      call h5pclose_f(crp_id, error)
      call h5dclose_f(dset_id, error)
      call h5fclose_f(file_id, error)
      !Close FORTRAN predefined datatypes
      call h5close_f(error)
    end subroutine
    subroutine GetData2D(file, dsetname, databuf)
      use HDF5
      character(len=*), intent(in) :: file, dsetname
      real, allocatable, intent(inout) :: databuf(:,:)

      integer(HSIZE_T) :: dimsr(2), maxdimsr(2)
      integer(HID_T) :: file_id, dspace_id, dset_id, crp_id, mem_id
      integer :: error, rankr

      !Init HDF
      call h5open_f(error)
      !Open the file.
      call h5fopen_f (file, H5F_ACC_RDONLY_F, file_id, error)
      !Open the  dataset.
      call h5dopen_f(file_id, dsetname, dset_id, error)
      !Get dataset's dataspace handle.
      call h5dget_space_f(dset_id, dspace_id, error)
      !Get dataspace's rank.
      call h5sget_simple_extent_ndims_f(dspace_id, rankr, error)
      !Get dataspace's dimensions.
      call h5sget_simple_extent_dims_f(dspace_id, dimsr, maxdimsr, error)
      !Get creation property list.
      CALL h5dget_create_plist_f(dset_id, crp_id, error)
      !create memory dataspace
      call h5screate_simple_f(rankr, dimsr, mem_id, error)
      !Read data
      allocate(databuf(dimsr(1),dimsr(2)))
      call H5dread_f(dset_id, H5T_NATIVE_DOUBLE, databuf, dimsr, &
           error, mem_id, dspace_id)
      !Close the objects that were opened.
      call h5sclose_f(dspace_id, error)
      call h5sclose_f(mem_id, error)
      call h5pclose_f(crp_id, error)
      call h5dclose_f(dset_id, error)
      call h5fclose_f(file_id, error)
      !Close FORTRAN predefined datatypes
      call h5close_f(error)
    end subroutine
end module
