program vfmimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar

  use srreheat, only : log_energy_reheat_ingev   
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use eosflow, only : eos_x, eos_norm_potential
  use eosflow, only : eos_norm_deriv_potential, eos_norm_deriv_second_potential
  use eosflow, only : eos_epsilon_one, eos_epsilon_two, eos_epsilon_three

  use vfmieos, only : vfmi_eos, vfmi_deriv_eos, vfmi_deriv_second_eos
  use vfmieos, only : vfmi_primitive_eos, vfmi_primitive_sqrteos

  use vfmisr, only : vfmi_epsilon_one, vfmi_epsilon_two,vfmi_epsilon_three
  use vfmisr, only : vfmi_norm_potential, vfmi_x_endinf, vfmi_numacc_betamax
  use vfmisr, only : vfmi_norm_deriv_potential, vfmi_norm_deriv_second_potential
  use vfmisr, only : vfmi_x_trajectory, vfmi_x_endinf
  use vfmireheat, only : vfmi_lnrhoreh_max, vfmi_x_star
  use vfmireheat, only : vfmi_x_rreh, vfmi_x_rrad

  use infinout, only : delete_file, livewrite
  
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh

  
  implicit none
  
  logical, parameter :: testParametric = .true.

  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts, n, nbeta

  integer, parameter :: nalpha = 6
  real(kp), dimension(nalpha), parameter :: alphavalues = (/0.5_kp,1.0_kp,1.5_kp,2.0_kp,2.5_kp,3._kp/)

  real(kp) :: bfold, bfoldMax, bfoldMin
  real(kp) :: wp1, dwp1, d2wp1, pwp1, psqrtwp1
  real(kp) :: x, xvfmi, V, dV, d2V

  real(kp) :: alpha, beta, w, bfoldstar
  real(kp) :: alphamin, alphamax, betamin, betamax
  real(kp) :: lnRhoReh, xstar, eps1, eps2, eps3, ns, r

  real(kp) :: lnRhoRehMin, lnRhoRehMax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar

  if (testParametric) then

     alpha = 1.9_kp
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
        
        psqrtwp1 = vfmi_primitive_sqrteos(bfold,alpha,beta)
        pwp1 = vfmi_primitive_eos(bfold,alpha,beta)
        wp1 = vfmi_eos(bfold,alpha,beta)
        dwp1 = vfmi_deriv_eos(bfold,alpha,beta)
        d2wp1 = vfmi_deriv_second_eos(bfold,alpha,beta)

        x = eos_x(psqrtwp1)
        xend = vfmi_x_endinf(alpha,beta)
        xvfmi = vfmi_x_trajectory(bfold,xend,alpha,beta)
        call livewrite('parametric_field.dat',bfold,x,xvfmi)

        V = eos_norm_potential(pwp1,wp1)
        dV = eos_norm_deriv_potential(pwp1,wp1,dwp1)
        d2V = eos_norm_deriv_second_potential(pwp1,wp1,dwp1,d2wp1)

        call livewrite('parametric_potential.dat',x,V,dV,d2V)

        V = vfmi_norm_potential(x,alpha,beta)
        dV = vfmi_norm_deriv_potential(x,alpha,beta)
        d2V = vfmi_norm_deriv_second_potential(x,alpha,beta)

        call livewrite('potential.dat',x,V,dV,d2V)

        eps1 = eos_epsilon_one(wp1)
        eps2 = eos_epsilon_two(wp1,dwp1)
        eps3 = eos_epsilon_three(wp1,dwp1,d2wp1)

        call livewrite('parametric_slowroll.dat',x,eps1,eps2,eps3)

        eps1 = vfmi_epsilon_one(x,alpha,beta)
        eps2 = vfmi_epsilon_two(x,alpha,beta)
        eps3 = vfmi_epsilon_three(x,alpha,beta)

        call livewrite('slowroll.dat',x,eps1,eps2,eps3)

     enddo

  end if

  npts = 20
  
  call aspicwrite_header('vfmi',labeps12,labnsr,labbfoldreh,(/'beta ','alpha'/))
  
  call delete_file('vfmi_predic.dat')
  call delete_file('vfmi_nsr.dat')

  do j=1,nalpha

     alpha = alphavalues(j)

     print *,'alpha= ',alpha

     nbeta = 50
     betamin=0.001
     betamax = min(100.,vfmi_numacc_betamax(300._kp,alpha)) 
     print *,'betamax=',betamax


     do k=1,nbeta
        beta = exp(log(betamin) + (log(betamax)-log(betamin))*real(k-1,kp)/real(nbeta-1,kp))

        lnRhoRehMin = lnRhoNuc
        xEnd = vfmi_x_endinf(alpha,beta)
        lnRhoRehMax = vfmi_lnrhoreh_max(alpha,beta,xend,Pstar)

        print *,'alpha= beta= lnRhoRehMin= lnRhoRehMax= ',alpha,beta,lnRhoRehMin,lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = vfmi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh= ',lnRhoReh, 'xstar= ', xstar, 'bfoldstar= ',bfoldstar

           eps1 = vfmi_epsilon_one(xstar,alpha,beta)
           eps2 = vfmi_epsilon_two(xstar,alpha,beta)
           eps3 = vfmi_epsilon_three(xstar,alpha,beta)

           logErehGev = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp-2._kp*eps1 - eps2
           r = 16._kp*eps1

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/beta,alpha/))
           
           !for sparing size
           if (abs(ns-1).gt.0.15) cycle
           if (r.lt.1e-10) cycle

           call livewrite('vfmi_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)
           call livewrite('vfmi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        enddo

     enddo

  enddo

  call aspicwrite_end()


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  n=2._kp
  alpha=0.5_kp
  beta = 0.1_kp
  xEnd = vfmi_x_endinf(alpha,beta)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = vfmi_x_rrad(alpha,beta,xend,lnRrad,Pstar,bfoldstar)
     print *
     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = vfmi_epsilon_one(xstar,alpha,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar

     eps1end =  vfmi_epsilon_one(xEnd,alpha,beta)
     VendOverVstar = vfmi_norm_potential(xEnd,alpha,beta)/vfmi_norm_potential(xstar,alpha,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = vfmi_x_rreh(alpha,beta,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = vfmi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo




end program vfmimain
