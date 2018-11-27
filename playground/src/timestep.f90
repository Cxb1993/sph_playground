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
  cndiff, cnhydro, difcond, mhdmuzero
integer, save ::&
  initialised=0,&
  n, eqSet(eqs_total)

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

  if (initialised == 0) call init()

  dt = huge(store(1,1))

  if (eqSet(eqs_diff)==1) then
    dt1 = cnhydro * minval(store(es_h,1:n)) / maxval(store(es_c,1:n))
    dt = min(dt, dt1)
  end if

  if (eqSet(eqs_magneto)==1) then
    dt1 = .1*minval(store(es_h,1:n))&
          /maxval(store(es_c,1:n))&
          /maxval(store(es_bx:es_bz,1:n))
    dt = min(dt, dt1)
  end if

  if (eqSet(eqs_diff)==1) then
    dt1 = cndiff * minval(store(es_den,1:n)) &
            * minval(store(es_c,1:n)) &
            * minval(store(es_h,1:n))**2 &
            / maxval(store(es_kappa,1:n))
    dt = min(dt, dt1)
  end if

  if (eqSet(eqs_fld)==1) then
    dt = min(dt, 1e-2)
  end if

  dtfinal = dt
end subroutine calc

end module
