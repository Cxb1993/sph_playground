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
      integer :: udump, umap
      real :: statetmp(ec_total)
      character(len=50) :: rfn, kfn

      call getState(statetmp)
      ! print*, statetmp
      call system_clock(start)

      call getresultfile(rfn)
      call getkerninflfilename(kfn)
      open(newunit=udump, file="output/newfulldump", status='replace', form='unformatted')
      open(newunit=umap, file="output/newdumpmap", status='replace', form='formatted')
      write(udump) t
      write(umap,*) 't 1 1 f8'
      write(udump) rfn, kfn
      write(umap,*) 'resultfile 1 1 50c'
      write(umap,*) 'influencefile 1 1 50c'
      write(udump) statetmp
      write(umap,*) 'state', ec_total, '1 f8'
      write(udump) store(1:es_total,1:int(statetmp(ec_realpn))+int(statetmp(ec_fixedpn)))
      write(umap,*) 'store', es_total, int(statetmp(ec_realpn))+int(statetmp(ec_fixedpn)), 'f8'
      close(udump)
      close(umap)
      call rename("output/newfulldump", "output/fulldump")
      call rename("output/newdumpmap", "output/dumpmap")

      call system_clock(finish)
      call addTime('dumper', finish - start)
    end subroutine dump

    subroutine restore(name, store, t)
      real, allocatable, intent(inout) :: store(:,:)
      real, intent(inout) :: t
      character(len=*), intent(in) :: name

      character(len=50) :: rfn, kfn
      integer :: u
      real :: statetmp(ec_total)

      call system_clock(start)
      print*, "#  # ---   ---   ---   ---   ---   ---   ---   ---   ---"
      open(newunit=u, file=name, action='read', form='unformatted')
      read(u) t
      write(*,blockFormatFlt) " #  #", " restored from "// name //" dump t = ", t
      read(u) rfn, kfn
      call setresultfile(rfn)
      call setkerninflfilename(kfn)
      read(u) statetmp
      call setState(statetmp)
      allocate(store(es_total,int(statetmp(ec_realpn))+int(statetmp(ec_fixedpn))))
      read(u) store
      close(u)
      print*, "#  # ---   ---   ---   ---   ---   ---   ---   ---   ---"
      call system_clock(finish)
      call addTime('dumper', finish - start)
    end subroutine restore

    subroutine clean()
      integer :: nu, snu
      open(newunit=nu, iostat=snu, file="output/newfulldump", status='old')
      ! open(newunit=du, iostat=sdu, file="output/fulldump", status='old')
      if (snu == 0) close(nu, status='delete')
      ! if (sdu == 0) close(du, status='delete')
      call rename("output/fulldump", "output/finaldump")
    end subroutine clean
end module
