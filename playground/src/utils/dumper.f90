module dumper
  use const
  use timing, only: addTime
  use state,  only: getState,&
                    setState,&
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
      real, allocatable, intent(in) :: store(:,:)
      real, intent(in) :: t
      integer :: u, rn, fn
      real :: statetmp(ec_total)
      character(len=50) :: rfn, kfn

      call getState(statetmp)
      ! print*, statetmp
      call system_clock(start)

      call getresultfile(rfn)
      call getkerninflfilename(kfn)
      open(newunit=u, file="output/newfulldump", status='replace')
      write(u,*) t
      write(u,*) rfn, kfn
      write(u,*) statetmp
      write(u,*) store(1:es_total,1:int(statetmp(ec_realpn))+int(statetmp(ec_fixedpn)))
      close(u)
      call rename("output/newfulldump", "output/fulldump")

      call system_clock(finish)
      call addTime('dumper', finish - start)
    end subroutine dump

    subroutine restore(store, t)
      real, allocatable, intent(inout) :: store(:,:)
      real, intent(inout) :: t
      character(len=50) :: rfn, kfn

      integer :: u
      real :: statetmp(ec_total)

      call system_clock(start)
      print*, "#  # ---   ---   ---   ---   ---   ---   ---   ---   ---"
      open(newunit=u, file="output/fulldump", action='read')
      read(u,*) t
      write(*,blockFormatFlt) " #  #", "restored from last dump t =", t
      read(u,*) rfn, kfn
      call setresultfile(rfn)
      call setkerninflfilename(kfn)
      read(u,*) statetmp
      call setState(statetmp)
      allocate(store(es_total,int(statetmp(ec_realpn))+int(statetmp(ec_fixedpn))))
      read(u,*) store
      close(u)
      print*, "#  # ---   ---   ---   ---   ---   ---   ---   ---   ---"
      call system_clock(finish)
      call addTime('dumper', finish - start)
    end subroutine restore

    subroutine clean()
      integer :: nu, du, snu, sdu
      open(newunit=nu, iostat=snu, file="output/newfulldump", status='old')
      open(newunit=du, iostat=sdu, file="output/fulldump", status='old')
      if (snu == 0) close(nu, status='delete')
      if (sdu == 0) close(du, status='delete')
    end subroutine clean
end module
