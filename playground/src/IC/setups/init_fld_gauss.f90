store(es_h,i) = hfac * sp
store(es_m,i) = (sp**dim) * rho1
store(es_den,i) = rho1
store(es_p,i) = prs1

store(es_t,i)  = 1./((2.*pi*(1./3.)**2)**(-dim/2.))*&
  exp(-0.5*(ra(1)*ra(1) + ra(2)*ra(2) + ra(3)*ra(3))/(1./3.)**2)
! store(es_u,i)  = store(es_t,i)
store(es_u,i) = prs1/(gamma - 1.)/rho1
store(es_kappa,i) = .1
