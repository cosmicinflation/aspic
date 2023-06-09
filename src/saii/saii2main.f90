program saii2main
  use infprec, only : kp, pi
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use saii2sr, only : saii2_epsilon_one, saii2_epsilon_two, saii2_epsilon_three
  use saii2reheat, only : saii2_lnrhoreh_max, saii2_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use saii2sr, only : saii2_norm_potential, saii2_x_endinf
  use saii2sr, only : saii2_norm_deriv_potential, saii2_norm_deriv_second_potential
  use saii2sr, only : saii2_numacc_mumin, saii2_numacc_efoldmax
  use saii2reheat, only : saii2_x_rreh, saii2_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use saiicommon, only : saii_x_potmax, saii_x_potzero, saii_x_epsoneunity
  
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  

  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,n,k
  integer :: npts = 15


  integer, parameter :: Na = 3
  real(kp) :: alphamin
  real(kp) :: alphamax

  integer, parameter :: Nmu=100
  real(kp) :: mumin = 0.1_kp
  real(kp) :: mumax = 100._kp
  
  real(kp) :: xmin = 0
  real(kp) :: xmax = 7
  
  real(kp) :: alpha,mu,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax, efoldmax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: V1,x,V,xvmax,xvzerop, xvzero, xvzerom

  real(kp), dimension(2) :: xeps
  
  Pstar = powerAmpScalar

  call delete_file('saii2_xV.dat')
  call delete_file('saii_xeps.dat')

  
  mu=1._kp

  n=250

  alphamin = -1.2_kp
  alphamax = 1.2_kp
  
  do i=1,n
     alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(n-1,kp)        

     xvmax = saii_x_potmax(alpha,mu)
     xvzero = saii_x_potzero(alpha,mu)

     if (alpha.le.-0.5_kp) then
        xvzerom = xvzero
        xvzerop = 2._kp*pi
     elseif (alpha.le.0._kp) then
        xvzerom = 0._kp
        xvzerop = 2._kp*pi
     else
        xvzerom = 0._kp
        xvzerop = xvzero
     endif
     
     call livewrite('saii2_xV.dat',alpha,xvmax,xvzerop,xvzerom)

     xeps = saii_x_epsoneunity(alpha,mu=10._kp)
     
     call livewrite('saii_xeps.dat',alpha,xeps(1),xeps(2))

  enddo

  alphamin = -0.8_kp
  alphamax = 0.6_kp

  
  call delete_file('saii2_predic.dat')
  call delete_file('saii2_nsr.dat')


  call aspicwrite_header('saii2',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha'/))
  
  w=0._kp

  do k=1,Na
     alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)
    
     mumin = saii2_numacc_mumin(120._kp,alpha)
 
     
     do j=0,Nmu 
        
        mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

        efoldmax = saii2_numacc_efoldmax(alpha,mu)
        print *,'alpha, mu= efoldNumAccMax= ',alpha,mu,efoldmax

        lnRhoRehMin = lnRhoNuc

        xend = saii2_x_endinf(alpha,mu)
        
        lnRhoRehMax = saii2_lnrhoreh_max(alpha,mu,xend,Pstar)

        print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=npts,1,-1

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = saii2_x_star(alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

           eps1 = saii2_epsilon_one(xstar,alpha,mu)
           eps2 = saii2_epsilon_two(xstar,alpha,mu)
           eps3 = saii2_epsilon_three(xstar,alpha,mu)


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)


           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('saii2_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('saii2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/mu,alpha/))

        end do

     end do
  enddo

  call aspicwrite_end()
  
  write(*,*)
  write(*,*)'Testing Rrad/Rreh'


  
  lnRradmin=-42
  lnRradmax = 10
  alpha = 0.7
  mu = 10._kp*saii2_numacc_mumin(120._kp,alpha)
  xend = saii2_x_endinf(alpha,mu)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = saii2_x_rrad(alpha,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = saii2_epsilon_one(xstar,alpha,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  saii2_epsilon_one(xend,alpha,mu)
     VendOverVstar = saii2_norm_potential(xend,alpha,mu) &
          /saii2_norm_potential(xstar,alpha,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = saii2_x_rreh(alpha,mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = saii2_x_star(alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program saii2main
