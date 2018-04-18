!test the reheating derivation from slow-roll
program ccsi3main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ccsi3sr, only : ccsi3_epsilon_one, ccsi3_epsilon_two, ccsi3_epsilon_three
  use ccsi3reheat, only : ccsi3_lnrhoreh_max, ccsi3_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use ccsicommon, only : ccsi_xmax
  use ccsi3sr, only : ccsi3_norm_potential, ccsi3_x_endinf, ccsi3_alphamin
  use ccsi3reheat, only : ccsi3_x_rreh, ccsi3_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat


  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,n
  integer :: npts = 20

  integer :: Np=10

  real(kp) :: alphamin= -1d-2
  real(kp) :: alphamax= -1d-6
  real(kp) :: efold

  real(kp) :: xmin = -2
  real(kp) :: xmax
  
  
  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: x, V3
  
  Pstar = powerAmpScalar



  call delete_file('ccsi3_potential.dat')
  call delete_file('ccsi3_slowroll.dat')


  alpha = -1e-4
  n=250

  xmax = ccsi_xmax(alpha)
  
  do i=1,n
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)        
    

     V3 = ccsi3_norm_potential(x,alpha)
     call livewrite('ccsi3_potential.dat',x*sqrt(1.5_kp),V3)

        
     eps1 = ccsi3_epsilon_one(x,alpha)
     eps2 = ccsi3_epsilon_two(x,alpha)
     eps3 = ccsi3_epsilon_three(x,alpha)
     call livewrite('ccsi3_slowroll.dat',x*sqrt(1.5_kp),eps1,eps2,eps3)

  enddo

  









  
  call delete_file('ccsi3_predic.dat')
  call delete_file('ccsi3_nsr.dat')

  efold = 120._kp

  alphamin = ccsi3_alphamin(efold)

  print *,'alphamin',alphamin

  !  w = 1._kp/3._kp
  w=0._kp

  do j=0,Np 

     alpha=-10._kp**(log10(-alphamin)+(log10(alphamax/alphamin))*(real(j,kp)/Np))

     print *

     print *,'alpha= ',alpha    

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = ccsi3_lnrhoreh_max(alpha,Pstar)

     print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = ccsi3_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)


        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = ccsi3_epsilon_one(xstar,alpha)
        eps2 = ccsi3_epsilon_two(xstar,alpha)
        eps3 = ccsi3_epsilon_three(xstar,alpha)


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)


        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('ccsi3_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('ccsi3_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 0.5*alphamin
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = ccsi3_x_rrad(alpha,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = ccsi3_epsilon_one(xstar,alpha)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = ccsi3_x_endinf(alpha)
     eps1end =  ccsi3_epsilon_one(xend,alpha)
     VendOverVstar = ccsi3_norm_potential(xend,alpha)/ccsi3_norm_potential(xstar,alpha)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = ccsi3_x_rreh(alpha,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = ccsi3_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program ccsi3main
