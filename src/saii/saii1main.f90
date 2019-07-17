program saii1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use saii1sr, only : saii1_epsilon_one, saii1_epsilon_two, saii1_epsilon_three
  use saii1reheat, only : saii1_lnrhoreh_max, saii1_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use saii1sr, only : saii1_norm_potential, saii1_x_endinf
  use saii1sr, only : saii1_norm_deriv_potential, saii1_norm_deriv_second_potential
  use saii1sr, only : saii1_numacc_mumin, saii1_numacc_efoldmax
  use saii1reheat, only : saii1_x_rreh, saii1_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use saiicommon, only : saii_x_potmax
  
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  

  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,n,k
  integer :: npts = 15


  integer, parameter :: Na = 3
  real(kp), parameter :: alphamin=-0.8_kp
  real(kp), parameter :: alphamax=0.6_kp

  integer, parameter :: Nmu=99
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

  real(kp) :: V1,x,V,V2


  
  Pstar = powerAmpScalar


  call delete_file('saii1_potential.dat')
  call delete_file('saii1_slowroll.dat')
  
  alpha = -0.2_kp
  mu=1._kp

  n=250

  
  do i=1,n
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)        

     V = saii1_norm_potential(x,alpha=0._kp,mu=1._kp)
     V1 = saii1_norm_potential(x,alpha=-0.2_kp,mu=1._kp)
     V2 = saii1_norm_potential(x,alpha=-0.8_kp,mu=1._kp)

     call livewrite('saii1_potential.dat',x,V1,V,V2)
     
     eps1 = saii1_epsilon_one(x,alpha,mu=10._kp)
     eps2 = saii1_epsilon_two(x,alpha,mu=10._kp)
     eps3 = saii1_epsilon_three(x,alpha,mu=10._kp)
          
     call livewrite('saii1_slowroll.dat',x,eps1,eps2,eps3)

     
  enddo

 
  call delete_file('saii1_predic.dat')
  call delete_file('saii1_nsr.dat')


  call aspicwrite_header('saii1',labeps12,labnsr,labbfoldreh,(/'mu   ','alpha'/))
  
  w=0._kp

  do k=1,Na
     alpha = alphamin +  real(k-1,kp)*(alphamax - alphamin)/real(Na-1,kp)
     
     mumin = saii1_numacc_mumin(120._kp,alpha)

     
     do j=0,Nmu 
        
        mu=10._kp**(log10(mumin)+(log10(mumax/mumin))*(real(j,kp)/real(Nmu,kp)))

        efoldmax = saii1_numacc_efoldmax(alpha,mu)
        print *,'alpha, mu= efoldNumAccMax= ',alpha,mu,efoldmax

        lnRhoRehMin = lnRhoNuc

        xend = saii1_x_endinf(alpha,mu)
        
        lnRhoRehMax = saii1_lnrhoreh_max(alpha,mu,xend,Pstar)

        print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=npts,1,-1

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = saii1_x_star(alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

           eps1 = saii1_epsilon_one(xstar,alpha,mu)
           eps2 = saii1_epsilon_two(xstar,alpha,mu)
           eps3 = saii1_epsilon_three(xstar,alpha,mu)


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)


           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('saii1_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('saii1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

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
  mu = 10._kp*saii1_numacc_mumin(120._kp,alpha)
  xend = saii1_x_endinf(alpha,mu)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = saii1_x_rrad(alpha,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = saii1_epsilon_one(xstar,alpha,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  saii1_epsilon_one(xend,alpha,mu)
     VendOverVstar = saii1_norm_potential(xend,alpha,mu) &
          /saii1_norm_potential(xstar,alpha,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = saii1_x_rreh(alpha,mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = saii1_x_star(alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program saii1main
