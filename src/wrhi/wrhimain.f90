!test the reheating derivation from slow-roll
program wrhimwrhin
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use wrhisr, only : wrhi_epsilon_one, wrhi_epsilon_two, wrhi_epsilon_three,wrhi_x_endinf
  use wrhireheat, only : wrhi_lnrhoreh_max, wrhi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use wrhisr, only : wrhi_norm_potential, wrhi_x_endinf
  use wrhireheat, only : wrhi_x_rreh, wrhi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  integer :: Nphi0
  real(kp) :: phi0min=10._kp**(-3.)
  real(kp) :: phi0max=10._kp**(3.)

  real(kp) :: phi0,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Nphi0=10

  Pstar = powerAmpScalar

  call delete_file('wrhi_predic.dat')
  call delete_file('wrhi_nsr.dat')

!  w = 1._kp/3._kp
  w=0._kp

 do j=0,Nphi0

 phi0=phi0min+(phi0max-phi0min)*(real(j,kp)/real(Nphi0,kp)) !arithmetic step
 phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(Nphi0,kp)) !logarithmic step

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = wrhi_lnrhoreh_max(phi0,Pstar)

  print *,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

	xstar = wrhi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar
 
       eps1 = wrhi_epsilon_one(xstar,phi0)
       eps2 = wrhi_epsilon_two(xstar,phi0)
       eps3 = wrhi_epsilon_three(xstar,phi0)

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)

       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('wrhi_predic.dat',phi0,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('wrhi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          Testing Rrad/Rreh           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  phi0 = 1._kp
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = wrhi_x_rrad(phi0,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = wrhi_epsilon_one(xstar,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = wrhi_x_endinf(phi0)
     eps1end =  wrhi_epsilon_one(xend,phi0)
     VendOverVstar = wrhi_norm_potential(xend,phi0)/wrhi_norm_potential(xstar,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = wrhi_x_rreh(phi0,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = wrhi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program wrhimwrhin
