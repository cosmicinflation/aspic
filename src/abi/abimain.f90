program abimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar

  use srreheat, only : log_energy_reheat_ingev   
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use eosflow, only : eos_x, eos_norm_potential
  use eosflow, only : eos_deriv_norm_potential, eos_deriv_second_norm_potential
  use eosflow, only : eos_epsilon_one, eos_epsilon_two, eps_epsilon_three

  use abieos, only : abi_eos, abi_deriv_eos, abi_deriv_second_eos
  use abieos, only : abi_eos_primitive, abi_sqrteos_primitive

  use abisr, only : abi_epsilon_one, abi_epsilon_two,abi_epsilon_three
  use abisr, only : abi_norm_potential, abi_x_endinf
  use abireheat, only : abi_lnrhoreh_max, abi_x_star
  use abireheat, only : abi_x_rreh, abi_x_rrad

  use infinout, only : delete_file, livewrite
  
  logical, parameter :: testParametric = .true.

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts, n, nalpha, nbeta

  real(kp) :: bfold, bfoldMax, bfoldMin
  real(kp) :: wp1, dwp1, d2wp1, pwp1, psqrtwp1
  real(kp) :: V, dV, d2V

  real(kp) :: alpha, beta, w, bfoldstar
  real(kp) :: alphamin, alphamax, betamin, betamax
  real(kp) :: lnRhoReh, xstar, eps1, eps2, eps3, ns, r

  real(kp) :: lnRhoRehMin, lnRhoRehMax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar

  if (testParametric) then

     alpha = 0.5_kp
     beta = 0.8_kp

     bfoldMin = -140._kp
     n=100

     call delete_file('parametric_potential.dat')
     call delete_file('parametric_slowroll.dat')
     call delete_file('parametric_field.dat')

     do i=1,n
        bfold = bfoldMin - real(i-1,kp)*bfoldMin/real(n-1,kp)
        
        pwp1 = abi_primitive_eos(bfold,alpha,beta)
        wp1 = abi_eos(bfold,alpha,beta)
        dwp1 = abi_deriv_eos(bfold,alpha,beta)
        d2wp1 = abi_deriv_second_eos(bfold,alpha,beta)

        V = eos_norm_potential(pwp1,wp1)
        dV = eos_norm_deriv_potential(pwp1,wp1,dwp1)
        d2V = eos_norm_deriv_second_potential(pwp1,wp1,dwp1,d2wp1)

        call livewrite('parametric_potential.dat',bfold,V,dV,d2V,d3V)

        
        eps1 = eos_epsilon_one(wp1)
        eps2 = eos_epsilon_two(wp1,dwp1)
        eps3 = eos_epsilon_three(wp1,dwp1,d2wp1)

        call livewrite('parametric_slowroll.dat',bfold,eps1,eps2,eps3)

        psqrtwp1 = abi_sqrteos_primitive(bfold,alpha,beta)

        x = eos_x(psqrtwp1)        
        call livewrite('parametric_field.dat',bdold,x)

     enddo

  end if




end program abimain
