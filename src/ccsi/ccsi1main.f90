!test the reheating derivation from slow-roll
program ccsi1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ccsi1sr, only : ccsi1_epsilon_one, ccsi1_epsilon_two, ccsi1_epsilon_three
  use ccsi1reheat, only : ccsi1_lnrhoreh_max, ccsi1_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use ccsi1sr, only : ccsi1_norm_potential, ccsi1_x_endinf
  use ccsi1sr, only : ccsi1_norm_deriv_potential, ccsi1_norm_deriv_second_potential
  use ccsi1reheat, only : ccsi1_x_rreh, ccsi1_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  

  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,n
  integer :: npts = 20

  integer :: Np=15
  real(kp) :: alphamin=1d-6
  real(kp) :: alphamax=1d-2

  real(kp) :: xmin = -2
  real(kp) :: xmax = 20
  
  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: V1,x,V
  
  Pstar = powerAmpScalar



  call delete_file('ccsi1_potential.dat')
  call delete_file('ccsi1_slowroll.dat')

  
  alpha = 1e-4
  n=250


  do i=1,n
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)        

     V = ccsi1_norm_potential(x,alpha=0._kp)
     V1 = ccsi1_norm_potential(x,alpha)
        

     call livewrite('ccsi1_potential.dat',x*sqrt(1.5_kp),V1,V)

     eps1 = ccsi1_epsilon_one(x,alpha)
     eps2 = ccsi1_epsilon_two(x,alpha)
     eps3 = ccsi1_epsilon_three(x,alpha)

     call livewrite('ccsi1_slowroll.dat',x*sqrt(1.5_kp),eps1,eps2,eps3)

  enddo

  


  call delete_file('ccsi1_predic.dat')
  call delete_file('ccsi1_nsr.dat')


  call aspicwrite_header('ccsi1',labeps12,labnsr,labbfoldreh,(/'alpha'/))
  
  !  w = 1._kp/3._kp
  w=0._kp

  do j=0,Np 

     alpha=10._kp**(log10(alphamin)+(log10(alphamax/alphamin))*(real(j,kp)/Np))

     print *

     print *,'alpha= ',alpha    

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = ccsi1_lnrhoreh_max(alpha,Pstar)

     print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = ccsi1_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)


        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = ccsi1_epsilon_one(xstar,alpha)
        eps2 = ccsi1_epsilon_two(xstar,alpha)
        eps3 = ccsi1_epsilon_three(xstar,alpha)


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)


        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('ccsi1_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('ccsi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha/))
        
     end do

  end do

  call aspicwrite_end()
  
  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 0.001
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = ccsi1_x_rrad(alpha,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = ccsi1_epsilon_one(xstar,alpha)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = ccsi1_x_endinf(alpha)
     eps1end =  ccsi1_epsilon_one(xend,alpha)
     VendOverVstar = ccsi1_norm_potential(xend,alpha)/ccsi1_norm_potential(xstar,alpha)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = ccsi1_x_rreh(alpha,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = ccsi1_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program ccsi1main
