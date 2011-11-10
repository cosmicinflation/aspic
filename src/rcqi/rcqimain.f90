!test the reheating derivation from slow-roll
program rcqimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rcqisr, only : rcqi_epsilon_one, rcqi_epsilon_two, rcqi_epsilon_three
  use rcqireheat, only : rcqi_lnrhoend, rcqi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i
  integer :: npts = 20

  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  Pstar = powerAmpScalar

!  call delete_file('rcqi_predic.dat')
!  call delete_file('rcqi_nsr.dat')

  alpha = 10._kp**(-0.485) 
!  w = 1._kp/3._kp
  w=0._kp
 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = rcqi_lnrhoend(alpha,Pstar)

  print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = rcqi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

       eps1 = rcqi_epsilon_one(xstar,alpha)
       eps2 = rcqi_epsilon_two(xstar,alpha)
       eps3 = rcqi_epsilon_three(xstar,alpha)

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('rcqi_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('rcqi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do



end program rcqimain
