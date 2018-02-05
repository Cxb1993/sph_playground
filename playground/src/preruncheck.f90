module preruncheck
  use const
  use state,  only: getdiffisotropic, &
                    getdiffconductivity, &
                    getmhdmagneticpressure,&
                    getPartNumber

  implicit none

  public checkVarsReady

contains
  subroutine checkVarsReady(s_tt, store)
    integer, intent(in) :: &
      s_tt
    real, allocatable, intent(in) :: &
      store(:,:)
    integer :: &
      difiso, i, rpn, fpn
    real :: &
      difcond, mhdmprs

    call getdiffisotropic(difiso)
    call getdiffconductivity(difcond)
    call getmhdmagneticpressure(mhdmprs)
    call getPartNumber(rpn, fpn)

    if ((mhdmprs <= 0.).and.((s_tt==eeq_magnetohydro).or.(s_tt==eeq_magnetohydrodiffusion))) then
      print* , "# # <!> mhd_magneticconstant <= 0"
      stop
    end if
    if ((difiso > 0).and.(difcond <= 0)) then
      print* , "# # <!> diff_conductivity <= 0"
      stop
    end if

    do i = 1,rpn+fpn
      if (int(store(es_type,i)) /= ept_empty) then
        if (store(es_den,i) <= 0.) then
          print*, "# # <!> density <= 0"
          stop
        end if
        if (store(es_c,i) <= 0.) then
          print*, "# # <!> csound <= 0"
          stop
        end if
        if (store(es_h,i) <= 0.) then
          print*, "# # <!> smoothing length <= 0"
          stop
        end if
        if (store(es_p,i) <= 0.) then
          print*, "# # <!> pressure <= 0"
          stop
        end if
        if (store(es_u,i) <= 0.) then
          print*, "# # <!> internal energy <= 0"
          stop
        end if
      end if
    end do
  end subroutine
end module
