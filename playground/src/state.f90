module state

  use const

  implicit none

  public :: set_difftype, get_difftype, set_tasktype, get_tasktype, &
            setddwtype, getddwtype, getdim, setdim, sinitvar, &
            ginitvar, setAdvancedDensity, getAdvancedDensity,&
            setArtificialTerms, getArtificialTerms, &
            scoordsys, gcoordsys, sorigin, gorigin, &
            setdiffconductivity, getdiffconductivity, &
            setdiffisotropic, getdiffisotropic, &
            setmhdmagneticpressure, getmhdmagneticpressure

  private
  save
    integer :: &
      dim = 1, ttype, ddwtype, dtype, icvar=-1, adden = 1, &
      artts = 1, coordsys = 1, origin = 0, diff_isotropic = 0
    real :: &
      diff_conductivity = 0.,&
      mhd_magneticconstant = 1.

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
      ! equations
      character (len=*), intent(in) :: itt
      select case(itt)
      case('hydro')
        ttype = eeq_hydro
      case('magnetohydro')
        ttype = eeq_magnetohydro
      case('diffusion')
        ttype = eeq_diffusion
      case('mhd')
        ttype = eeq_magnetohydrodiffusion
      ! case('diff-laplace')
      !   ttype = 5
      ! case('diff-graddiv')
      !   ttype = 6
      ! case('chi-laplace')
      !   ttype = 7
      ! case('chi-graddiv')
      !   ttype = 8
      ! case('diff-artvisc')
      !   ttype = 10
      case default
        print *, 'Task type not set: ', itt
        stop
      end select
    end subroutine set_tasktype

    pure subroutine get_tasktype(ott)
      integer, intent(out) :: ott
      ott = ttype
    end subroutine get_tasktype

    subroutine setddwtype(itt)
      character (len=*), intent(in) :: itt
      select case(itt)
      case('n2w')
        ddwtype = esd_n2w
      case('fab')
        ddwtype = esd_fab
      case('2nw')
        ddwtype = esd_2nw
      case('fw')
        ddwtype = esd_fw
      case default
        ddwtype = esd_fab
        write(*,blockFormatStr2) ' # <?>', ' Default d2W/dx2: ', 'fab'
      end select
     !  call calc_params()
    end subroutine setddwtype

    pure subroutine getddwtype(ott)
      integer, intent(out) :: ott
      ott = ddwtype
    end subroutine getddwtype

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
     case('otvortex')
       icvar = ett_OTvortex
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
       write(*,blockFormatStr2) ' # <?>', ' Default artificial terms: ', 'yes'
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
       write(*,blockFormatStr2) ' # <?>', ' Default coordinate system: ', 'cartesian'
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

   subroutine setdiffconductivity(i)
     real, intent(in) :: i
     diff_conductivity = i
   end subroutine
   pure subroutine getdiffconductivity(o)
     real, intent(out) :: o
     o = diff_conductivity
   end subroutine

   subroutine setdiffisotropic(i)
     integer, intent(in) :: i
     diff_isotropic = i
   end subroutine
   pure subroutine getdiffisotropic(o)
     integer, intent(out) :: o
     o = diff_isotropic
   end subroutine

   subroutine setmhdmagneticpressure(i)
     real, intent(in) :: i
     mhd_magneticconstant = i
   end subroutine
   pure subroutine getmhdmagneticpressure(o)
     real, intent(out) :: o
     o = mhd_magneticconstant
   end subroutine
end module
