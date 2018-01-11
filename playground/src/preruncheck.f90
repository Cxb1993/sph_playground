module preruncheck
  use const
  use state,  only: getdiffisotropic, &
                    getdiffconductivity, &
                    getmhdmagneticpressure
  use BC,     only: realpartnumb

  implicit none

  public checkVarsReady

contains
  subroutine checkVarsReady(s_tt, store)
    integer, intent(in) :: &
      s_tt
    real, allocatable, intent(in) :: &
      store(:,:)
    integer :: &
      err = 0, difiso
    real :: &
      difcond, mhdmprs

    call getdiffisotropic(difiso)
    call getdiffconductivity(difcond)
    call getmhdmagneticpressure(mhdmprs)

    if (minval(store(es_den,1:realpartnumb)) <= 0.) then
      print*, "# # <!> density <= 0"
      err = 1
    end if
    if (minval(store(es_c,1:realpartnumb)) <= 0.) then
      print*, "# # <!> csound <= 0"
      err = 1
    end if
    if (minval(store(es_h,1:realpartnumb)) <= 0.) then
      print*, "# # <!> smoothing length <= 0"
      err = 1
    end if
    if (minval(store(es_p,1:realpartnumb)) <= 0.) then
      print*, "# # <!> pressure <= 0"
      err = 1
    end if
    if (minval(store(es_u,1:realpartnumb)) <= 0.) then
      print*, "# # <!> internal energy <= 0"
      err = 1
    end if
    if ((mhdmprs <= 0.).and.((s_tt==eeq_magnetohydro).or.(s_tt==eeq_magnetohydrodiffusion))) then
      print* , "# # <!> mhd_magneticconstant <= 0"
      err = 1
    end if
    if ((difiso > 0).and.(difcond <= 0)) then
      print* , "# # <!> diff_conductivity <= 0"
      err = 1
    end if
    if (err == 1) then
      print*, "# # (Ш_Ш) preruncheck: DOOOMED"
      stop
    end if
  end subroutine
end module
