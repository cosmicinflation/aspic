program abimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar

  use srreheat, only : log_energy_reheat_ingev   
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use eosflow, only : eos_absx, eos_norm_potential
  use eosflow, only : eos_norm_absderiv_potential, eos_norm_deriv_second_potential
  use eosflow, only : eos_epsilon_one, eos_epsilon_two, eos_epsilon_three

  use abieos, only : abi_eos, abi_deriv_eos, abi_deriv_second_eos
  use abieos, only : abi_eos_primitive, abi_sqrteos_primitive

  use abisr, only : abi_epsilon_one, abi_epsilon_two,abi_epsilon_three
  use abisr, only : abi_norm_potential, abi_x_endinf, abi_numacc_betamax
  use abisr, only : abi_norm_deriv_potential, abi_norm_deriv_second_potential
  use abireheat, only : abi_lnrhoreh_max, abi_x_star
  use abireheat, only : abi_x_rreh, abi_x_rrad

  use infinout, only : delete_file, livewrite

  implicit none
  
  logical, parameter :: testParametric = .true.

  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts, n, nbeta

  integer, parameter :: nalpha = 6
  real(kp), dimension(nalpha), parameter :: alphavalues = (/0.5_kp,1.0_kp,1.5_kp,2.0_kp,2.5_kp,3._kp/)

  real(kp) :: bfold, bfoldMax, bfoldMin
  real(kp) :: wp1, dwp1, d2wp1, pwp1, psqrtwp1
  real(kp) :: x, V, dV, d2V

  real(kp) :: alpha, beta, w, bfoldstar
  real(kp) :: alphamin, alphamax, betamin, betamax
  real(kp) :: lnRhoReh, xstar, eps1, eps2, eps3, ns, r

  real(kp) :: lnRhoRehMin, lnRhoRehMax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar

  if (testParametric) then

     alpha = 2.1_kp
     beta = 0.8_kp

     bfoldMin = -140._kp
     n=250

     call delete_file('parametric_potential.dat')
     call delete_file('parametric_slowroll.dat')
     call delete_file('parametric_field.dat')

     call delete_file('potential.dat')
     call delete_file('slowroll.dat')

     do i=1,n
        bfold = bfoldMin - real(i-1,kp)*bfoldMin/real(n-1,kp)        
        
        psqrtwp1 = abi_sqrteos_primitive(bfold,alpha,beta)
        pwp1 = abi_eos_primitive(bfold,alpha,beta)
        wp1 = abi_eos(bfold,alpha,beta)
        dwp1 = abi_deriv_eos(bfold,alpha,beta)
        d2wp1 = abi_deriv_second_eos(bfold,alpha,beta)

        x = eos_absx(psqrtwp1)
        call livewrite('parametric_field.dat',bfold,x)

        V = eos_norm_potential(pwp1,wp1)
        dV = eos_norm_absderiv_potential(pwp1,wp1,dwp1)
        d2V = eos_norm_deriv_second_potential(pwp1,wp1,dwp1,d2wp1)

        call livewrite('parametric_potential.dat',x,V,dV,d2V)

        V = abi_norm_potential(x,alpha,beta)
        dV = abi_norm_deriv_potential(x,alpha,beta)
        d2V = abi_norm_deriv_second_potential(x,alpha,beta)

        call livewrite('potential.dat',x,V,dV,d2V)

        eps1 = eos_epsilon_one(wp1)
        eps2 = eos_epsilon_two(wp1,dwp1)
        eps3 = eos_epsilon_three(wp1,dwp1,d2wp1)

        call livewrite('parametric_slowroll.dat',x,eps1,eps2,eps3)

        eps1 = abi_epsilon_one(x,alpha,beta)
        eps2 = abi_epsilon_two(x,alpha,beta)
        eps3 = abi_epsilon_three(x,alpha,beta)

        call livewrite('slowroll.dat',x,eps1,eps2,eps3)

     enddo

  end if


  call delete_file('abi_predic.dat')
  call delete_file('abi_nsr.dat')

  do j=1,nalpha

     alpha = alphavalues(j)

     print *,'alpha= ',alpha

     nbeta = 10
     betamin=0.001
     betamax = min(10.,abi_numacc_betamax(300._kp,alpha)) 
     print *,'betamax=',betamax



     npts = 20

     do k=1,nbeta
        beta = exp(log(betamin) + (log(betamax)-log(betamin))*real(k-1,kp)/real(nbeta-1,kp))

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = abi_lnrhoreh_max(alpha,beta,Pstar)

        print *,'alpha= beta= lnRhoRehMin= lnRhoRehMax= ',alpha,beta,lnRhoRehMin,lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = abi_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh= ',lnRhoReh, 'xstar= ', xstar, 'bfoldstar= ',bfoldstar

           eps1 = abi_epsilon_one(xstar,alpha,beta)
           eps2 = abi_epsilon_two(xstar,alpha,beta)
           eps3 = abi_epsilon_three(xstar,alpha,beta)

           logErehGev = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp-2._kp*eps1 - eps2
           r = 16._kp*eps1

           !for sparing size
           !           if (abs(ns-1).gt.0.15) cycle
           !           if (r.lt.1e-10) cycle

           call livewrite('abi_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)
           call livewrite('abi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        enddo

     enddo

  enddo

end program abimain
