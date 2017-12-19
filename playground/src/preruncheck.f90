module preruncheck

  public checkVarsReady

contains
  subroutine checkVarsReady(s_tt, mhd_magneticconstant, rho, csnd, sln, prs, intu)
    use const
    use BC, only: realpartnumb
    implicit none

    integer, intent(in) :: s_tt
    real, intent(in) :: mhd_magneticconstant
    real, allocatable, intent(in) :: rho(:), csnd(:), sln(:), prs(:), intu(:)
    integer :: err = 0
    
    if (minval(rho(1:realpartnumb)) <= 0.) then
      print*, "+++> density <= 0"
      err = 1
    end if
    if (minval(csnd(1:realpartnumb)) <= 0.) then
      print*, "+++> csound <= 0"
      err = 1
    end if
    if (minval(sln(1:realpartnumb)) <= 0.) then
      print*, "+++> smoothing length <= 0"
      err = 1
    end if
    if (minval(prs(1:realpartnumb)) <= 0.) then
      print*, "+++> pressure <= 0"
      err = 1
    end if
    if (minval(intu(1:realpartnumb)) <= 0.) then
      print*, "+++> internal energy <= 0"
      err = 1
    end if
    if ((mhd_magneticconstant <= 0.).and.((s_tt==eeq_magnetohydro).or.(s_tt==eeq_magnetohydrodiffusion))) then
      print* , "+++> mhd_magneticconstant <= 0"
      err = 1
    end if
    if (err == 1) then
      print*, "DOOOMED"
      stop
    end if
  end subroutine
end module
