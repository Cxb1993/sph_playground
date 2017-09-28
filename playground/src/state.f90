module state
  implicit none

  public :: set_difftype, get_difftype, set_tasktype, get_tasktype, &
            set_kerntype, get_kerntype, getdim, setdim, sinitvar, &
            ginitvar, setAdvancedDensity, getAdvancedDensity,&
            setArtificialTerms, getArtificialTerms
  private
  save
    integer :: dim = 1
    integer :: ttype, ktype, dtype, icvar=-1, adden = 1, artts = 1
  contains
    subroutine setdim(d)
      integer, intent(in) :: d
      dim = d
    end subroutine

    pure subroutine getdim(d)
      integer, intent(out) :: d
      d = dim
    end subroutine

    subroutine set_tasktype(itt)
      character (len=*), intent(in) :: itt
      select case(itt)
      case('hydroshock')
        ttype = 1
      case('infslb')
        ttype = 2
      case('heatconduction')
        ttype = 3
      case('pheva')
        ttype = 4
      case('diff-laplace')
        ttype = 5
      case('diff-graddiv')
        ttype = 6
      case('chi-laplace')
        ttype = 7
      case('chi-graddiv')
        ttype = 8
      case('soundwave')
        ttype = 9
      case('chi-artvisc')
        ttype = 10
      case default
        print *, 'Task type not set: ', itt
        stop
      end select
    end subroutine set_tasktype

    pure subroutine get_tasktype(ott)
      integer, intent(out) :: ott
      ott = ttype
    end subroutine get_tasktype

    subroutine set_kerntype(itt)
      character (len=*), intent(in) :: itt
      select case(itt)
      case('n2w')
        ktype = 1
      case('fab')
        ktype = 2
      case('2nw')
        ktype = 3
      case('fw')
        ktype = 4
      case default
        print *, 'Kernel type not set: ', itt
        stop
      end select
     !  call calc_params()
    end subroutine set_kerntype

    pure subroutine get_kerntype(ott)
      integer, intent(out) :: ott
      ott = ktype
    end subroutine get_kerntype

   subroutine set_difftype(idt)
     character (len=*), intent(in) :: idt
     select case(idt)
     case('diff')
       dtype = 1
     case('symm')
       dtype = 2
     case default
       print *, 'Differentiation type is not set: ', idt
       stop
     end select
   end subroutine set_difftype

   pure subroutine get_difftype(odt)
     integer, intent(out) :: odt
     odt = dtype
   end subroutine get_difftype

   subroutine sinitvar(itt)
     character (len=*), intent(in) :: itt
     select case(itt)
     case('isotropic-sinxsinysinz')
       icvar = 1
     case('anisotropic-sinxsinysinz')
       icvar = 2
     case default
       print *, 'There is no such initial variable setting : ', itt
       stop
     end select
   end subroutine

   pure subroutine ginitvar(ott)
     integer, intent(out) :: ott
     ott = icvar
   end subroutine

   subroutine setAdvancedDensity(iad)
     character (len=*), intent(in) :: iad
     select case(iad)
     case('yes')
       adden = 1
     case('no')
       adden = 0
     case default
       print *, 'There is no such case for advanced density : ', iad
       stop
     end select
   end subroutine

   pure subroutine getAdvancedDensity(oad)
     integer, intent(out) :: oad
     oad = adden
   end subroutine

   subroutine setArtificialTerms(iat)
     character (len=*), intent(in) :: iat
     select case(iat)
     case('yes')
       artts = 1
     case('no')
       artts = 0
     case default
       print *, 'There is no such case for advanced density : ', iat
       stop
     end select
   end subroutine

   pure subroutine getArtificialTerms(oat)
     integer, intent(out) :: oat
     oat = artts
   end subroutine
end module
