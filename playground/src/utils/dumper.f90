module dumper
  use const
  use timing, only: addTime
  use state,  only: setState,&
                    getresultfile,&
                    getkerninflfilename,&
                    setresultfile,&
                    setkerninflfilename

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

      allocate(state(ec_total))
      call GetData1D("output/dump_full.h5",&
        "state", ec_total, state)

      call setState(state)
      allocate(store(es_total,int(state(ec_realpn))+int(state(ec_fixedpn))))

      write(*,blockFormatFltSci) " #  #", " restored from t = ", state(ec_time)
      write(*,blockFormatFltSci) " #  #", "             dump #", state(ec_lastprint)

      call GetData2D("output/dump_full.h5",&
        "store", es_total, int(state(ec_realpn))+int(state(ec_fixedpn)), store)

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
      INTEGER(HSIZE_T) :: ddim(1)
      integer, intent(in) :: dim

      integer(HID_T) :: file_id, dspace_id, dset_id
      integer :: error

      call h5open_f(error)
      call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error)
      ddim(1) = dim
      call h5screate_simple_f(1, ddim, dspace_id, error)
      call h5dcreate_f(file_id, dspace, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, ddim, error)
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)
    end subroutine
    subroutine PutData2D(file, dspace, dim1, dim2, data)
      use HDF5
      character(len=*), intent(in) :: file, dspace
      real, allocatable, intent(in) :: data(:,:)
      INTEGER(HSIZE_T) :: ddim(2)
      integer, intent(in) :: dim1, dim2

      integer(HID_T) :: file_id, dspace_id, dset_id
      integer :: error

      call h5open_f(error)
      call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error)
      ddim(1) = dim1
      ddim(2) = dim2
      call h5screate_simple_f(2, ddim, dspace_id, error)
      call h5dcreate_f(file_id, dspace, H5T_NATIVE_DOUBLE, dspace_id, dset_id, error)
      call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, ddim, error)
      call h5dclose_f(dset_id, error)
      call h5sclose_f(dspace_id, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)
    end subroutine

    subroutine GetData1D(file, dspace, dim, data)
      use HDF5
      character(len=*), intent(in) :: file, dspace
      real, allocatable, intent(inout) :: data(:)
      integer, intent(in) :: dim

      INTEGER(HSIZE_T) :: ddim(1)
      integer(HID_T) :: file_id, dspace_id, dset_id
      integer :: error

      call h5open_f(error)
      call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error)
      ddim(1) = dim
      call h5dopen_f(file_id, dspace, dset_id, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, ddim, error)
      call h5dclose_f(dset_id, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)
    end subroutine
    subroutine GetData2D(file, dspace, dim1, dim2, data)
      use HDF5
      character(len=*), intent(in) :: file, dspace
      real, allocatable, intent(inout) :: data(:,:)
      integer, intent(in) :: dim1, dim2

      INTEGER(HSIZE_T) :: ddim(2)
      integer(HID_T) :: file_id, dspace_id, dset_id
      integer :: error

      call h5open_f(error)
      call h5fopen_f(file, H5F_ACC_RDWR_F, file_id, error)
      ddim(1) = dim1
      ddim(2) = dim2
      call h5dopen_f(file_id, dspace, dset_id, error)
      call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, data, ddim, error)
      call h5dclose_f(dset_id, error)
      call h5fclose_f(file_id, error)
      call h5close_f(error)
    end subroutine
end module
