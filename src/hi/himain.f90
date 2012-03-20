!test the reheating derivation from slow-roll
program himain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use hisr, only : hi_epsilon_one, hi_epsilon_two, hi_epsilon_three
  use hireheat, only : hi_lnrhoend, hi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i
  integer :: npts = 20

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer



  Pstar = powerAmpScalar

  call delete_file('hi_predic.dat')
  call delete_file('hi_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp


  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = hi_lnrhoend(Pstar)

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



end program himain
