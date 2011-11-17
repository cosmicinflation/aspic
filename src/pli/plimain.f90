!test the reheating derivation from slow-roll
program plimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use plisr, only : pli_epsilon_one, pli_epsilon_two, pli_epsilon_three
  use plireheat, only : pli_lnrhoend, pli_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i
  integer :: npts = 2

  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  Pstar = powerAmpScalar

  call delete_file('pli_predic.dat')
  call delete_file('pli_nsr.dat')

  alpha=sqrt(2._kp)*10**(-5./2.)
  do while (alpha.le.0.9_kp)

  alpha = alpha*1.1
!  w = 1._kp/3._kp
  w=0._kp
 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = pli_lnrhoend(alpha,Pstar)

  print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = pli_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

       eps1 = pli_epsilon_one(xstar,alpha)
       eps2 = pli_epsilon_two(xstar,alpha)
       eps3 = pli_epsilon_three(xstar,alpha)

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('pli_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('pli_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

  enddo



end program plimain
