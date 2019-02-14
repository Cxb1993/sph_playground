module timestep

use const
use kernel,           only: getcndiff,&
                            getcnhydro
use state,            only: getStateVal,&
                            getEqComponent

implicit none

public calc

private

real, save ::&
  cndiff, cnhydro, difcond, mhdmuzero, cnfld
integer, save ::&
  initialised=0,&
  n, eqSet(eqs_total)
real, parameter::&
  timestepcut = 1e-30

contains

subroutine init()
  real ::&
    nr
  call getcnhydro(cnhydro)
  call getcndiff(cndiff)
  call getStateVal(ec_dcondconst, difcond)
  call getStateVal(ec_muzero, mhdmuzero)
  call getStateVal(ec_realpn, nr)
  n = int(nr)
  call getEqComponent(eqSet)
  initialised = 1
end subroutine init

subroutine calc(store, dtfinal)
  real, allocatable, intent(in) :: &
    store(:,:)
  real, intent(out) ::&
    dtfinal
  real ::&
    dt, dt1

  real ::&
    lightspeed, nta(3), ka, da, ta, eddfact, fld_ra, lambdaa,&
    ha, maxeddfact, maxlambda, maxt, minden, minh, minkappa,&
    maxkappa, maxu, minksi, mint, ua, maxden, minu

  real::&
    alpha, cv, gamma, mu, Rg, sigmaB
  integer ::&
    i

  if (initialised == 0) call init()

  dt = huge(store(1,1))

  if (eqSet(eqs_hydro)==1) then
    dt1 = cnhydro * minval(store(es_h,1:n)) / maxval(store(es_c,1:n))
    ! print*, dt1
    dt = min(dt, dt1)
  end if

  if (eqSet(eqs_magneto)==1) then
    dt1 = .1*minval(store(es_h,1:n))&
          /maxval(store(es_c,1:n))&
          /maxval(store(es_bx:es_bz,1:n))
    ! print*, dt1
    dt = min(dt, dt1)
  end if

  if (eqSet(eqs_diff)==1) then
    dt1 = cndiff * minval(store(es_den,1:n)) &
            * minval(store(es_c,1:n)) &
            * minval(store(es_h,1:n))**2 &
            / maxval(store(es_kappa,1:n))
    ! print*, dt1
    dt = min(dt, dt1)
  end if

  if (eqSet(eqs_fld)==1) then
    lightspeed = 2.997924e10
    sigmaB = 5.6704e-5
    mu = 0.2
    Rg = 8.314e7
    gamma = 5./3.

    alpha = 4.0*sigmaB/lightspeed
    cv = (gamma-1)*mu/Rg

    minden = store(es_den,1)
    maxden = store(es_den,1)
    minkappa = store(es_kappa,1)
    maxkappa = store(es_kappa,1)
    minh = store(es_h,1)
    maxt = store(es_t,1)
    mint = store(es_t,1)
    maxu = store(es_u,1)
    minu = store(es_u,1)
    maxeddfact = 0.0
    maxlambda = 0.0
    minksi = maxt/maxden

    do i = 1, n
      nta = store(es_dtdx:es_dtdz,i)
      ka = store(es_kappa,i)
      da = store(es_den,i)
      ta = store(es_t,i)
      ha = store(es_h,i)
      ua = store(es_u,i)

      if (minden > da) minden = da
      if (maxden < da) maxden = da
      if (minkappa > ka) minkappa = ka
      if (maxkappa < ka) maxkappa = ka
      if (minh > ha) minh = ha
      if (maxu < ua) maxu = ua
      if (minu > ua) minu = ua
      if (maxt < ta) maxt = ta
      if (mint > ta) then
        mint = ta
        minksi = ta/da
      end if

      fld_ra = sqrt(dot_product(nta(:),nta(:)))/(ka*da*ta)
      lambdaa = (2. + fld_ra)/(6. + 3*fld_ra + fld_ra*fld_ra)
      eddfact = lambdaa + lambdaa*lambdaa*fld_ra*fld_ra

      if (maxlambda < lambdaa) maxlambda = lambdaa
      if (maxeddfact < eddfact) maxeddfact = eddfact
    end do

    cnfld = 0.3
    !
    ! dt1 = cnfld*minh*minh*minden*minkappa/lightspeed/maxlambda
    ! ! print*, dt1
    ! dt = min(dt, dt1)
    !
    dt1 = cnfld/maxeddfact/maxt
    ! print*, dt1
    dt = min(dt, dt1)

    ! dt1 = 0.3*minksi/(alpha*lightspeed*maxkappa*abs(maxt/alpha-(maxu/cv)**4))
    ! dt = min(dt, dt1)
    ! dt1 = cnfld*minksi/alpha/lightspeed/maxkappa/(maxt/alpha)
    ! print*, dt1
    ! if (dt1 > timestepcut) dt = min(dt, dt1)
    !
    ! dt1 = cnfld*minksi/alpha/lightspeed/maxkappa/(maxu/cv)**4
    ! print*, dt1
    ! if (dt1 > timestepcut) dt = min(dt, dt1)
    !
    !
    ! dt1 = cnfld*minu/alpha/lightspeed/maxkappa/(maxt/alpha)
    ! print*, dt1
    ! if (dt1 > timestepcut) dt = min(dt, dt1)
    !
    ! dt1 = cnfld*minu/alpha/lightspeed/maxkappa/(maxu/cv)**4
    ! print*, dt1
    ! if (dt1 > timestepcut) dt = min(dt, dt1)
    ! dt = 1e5
  end if
  ! print*, dt
  ! print*, '=============='
  ! read*

  dtfinal = dt
end subroutine calc

end module
