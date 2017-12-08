module state

  use const

  implicit none

  public :: set_difftype, get_difftype, set_tasktype, get_tasktype, &
            setkerntype, getkerntype, getdim, setdim, sinitvar, &
            ginitvar, setAdvancedDensity, getAdvancedDensity,&
            setArtificialTerms, getArtificialTerms, &
            setpartnum, getpartnum, scoordsys, gcoordsys, &
            sorigin, gorigin, switch_hc_conductivity, switch_hc_isotropic

  private
  save
    integer :: dim = 1, partnumber=-1
    integer :: ttype, ktype, dtype, icvar=-1, adden = 1, artts = 1, coordsys = 1, &
                origin = 0
    real :: &
      switch_hc_conductivity = 0.,&
      switch_hc_isotropic = -1.

  contains
    subroutine setdim(d)
      integer, intent(in) :: d
      dim = d
    end subroutine

    pure subroutine getdim(d)
      integer, intent(out) :: d
      d = dim
    end subroutine

    subroutine setpartnum(d)
      integer, intent(in) :: d
      partnumber = d
    end subroutine

    pure subroutine getpartnum(d)
      integer, intent(out) :: d
      d = partnumber
    end subroutine

    subroutine set_tasktype(itt)
      ! equations
      character (len=*), intent(in) :: itt
      select case(itt)
      case('hydro')
        ttype = 1
      case('magnetohydro')
        ttype = 2
      case('diffusion')
        ttype = 3
      case('hmd')
        ttype = 4
      case('diff-laplace')
        ttype = 5
      case('diff-graddiv')
        ttype = 6
      case('chi-laplace')
        ttype = 7
      case('chi-graddiv')
        ttype = 8
      case('diff-artvisc')
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

    subroutine setkerntype(itt)
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
        ktype = 2
        print *, '# # > Fab is set be default d2W/dx2 <'
      end select
     !  call calc_params()
    end subroutine setkerntype

    pure subroutine getkerntype(ott)
      integer, intent(out) :: ott
      ott = ktype
    end subroutine getkerntype

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
     ! initvar
     character (len=*), intent(in) :: itt
     select case(itt)
     case('sin3')
       icvar = ett_sin3
     case('mti')
       icvar = ett_mti
     case('shock12')
       icvar = ett_shock12
     case('pulse')
       icvar = ett_pulse
     case('ring')
       icvar = ett_ring
     case('soundwave')
       icvar = ett_soundwave
     case('hydroshock')
       icvar = ett_hydroshock
     case('alfvenwave')
       icvar = ett_alfvenwave
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
       artts = 1
       print *, '# # > Yes is set be default choise for artificial terms <'
     end select
   end subroutine

   pure subroutine getArtificialTerms(oat)
     integer, intent(out) :: oat
     oat = artts
   end subroutine

   subroutine scoordsys(ics)
     character (len=*), intent(in) :: ics
     select case(ics)
     case('cartesian')
       coordsys = 1
     case('cylindric')
       coordsys = 2
     case default
       coordsys = 1
       print *, '# # > Cartesian is set be default coordinate system <'
      end select
   end subroutine

   pure subroutine gcoordsys(ocs)
     integer, intent(out) :: ocs
     ocs = coordsys
   end subroutine

   subroutine sorigin(io)
     integer, intent(in) :: io
     origin = io
   end subroutine sorigin

   pure subroutine gorigin(oo)
     integer, intent(out) :: oo
     oo = origin
   end subroutine gorigin
end module
