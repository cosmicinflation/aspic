!test the reheating derivation from slow-roll
program himain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use hisr, only : hi_epsilon_one, hi_epsilon_two, hi_epsilon_three
  use hireheat, only : hi_lnrhoreh_max, hi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use hisr, only : hi_norm_potential, hi_x_endinf
  use hireheat, only : hi_x_rreh, hi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i
  integer :: npts = 20

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar

  call delete_file('hi_predic.dat')
  call delete_file('hi_nsr.dat')


  !  w = 1._kp/3._kp
  w=0._kp


  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = hi_lnrhoreh_max(Pstar)

  print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)



     xstar = hi_x_star(w,lnRhoReh,Pstar,bfoldstar)



     print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar


     eps1 = hi_epsilon_one(xstar)
     eps2 = hi_epsilon_two(xstar)
     eps3 = hi_epsilon_three(xstar)


     logErehGeV = log_energy_reheat_ingev(lnRhoReh)


     Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


     ns = 1._kp - 2._kp*eps1 - eps2
     r =16._kp*eps1

     call livewrite('hi_predic.dat',eps1,eps2,eps3,r,ns,Treh)

     call livewrite('hi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('hi_predic_summarized.dat')
  lnRhoReh = lnRhoRehMin
  xstar = hi_x_star(w,lnRhoReh,Pstar,bfoldstar)
  eps1 = hi_epsilon_one(xstar)
  eps2 = hi_epsilon_two(xstar)
  eps3 = hi_epsilon_three(xstar)
  ns = 1._kp - 2._kp*eps1 - eps2
  r =16._kp*eps1
  call livewrite('hi_predic_summarized.dat',eps1,eps2,eps3,r,ns)
  lnRhoReh = lnRhoRehMax
  xstar = hi_x_star(w,lnRhoReh,Pstar,bfoldstar)
  eps1 = hi_epsilon_one(xstar)
  eps2 = hi_epsilon_two(xstar)
  eps3 = hi_epsilon_three(xstar)
  ns = 1._kp - 2._kp*eps1 - eps2
  r =16._kp*eps1
  call livewrite('hi_predic_summarized.dat',eps1,eps2,eps3,r,ns)

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = hi_x_rrad(lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = hi_epsilon_one(xstar)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = hi_x_endinf()
     eps1end =  hi_epsilon_one(xend)
     VendOverVstar = hi_norm_potential(xend)/hi_norm_potential(xstar)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = hi_x_rreh(lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = hi_x_star(w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program himain
