!test the reheating derivation from slow-roll
program aimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use aisr, only : ai_epsilon_one, ai_epsilon_two, ai_epsilon_three,ai_x_endinf
  use aireheat, only : ai_lnrhoend, ai_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh, mu, mumin, mumax

  integer :: i,j
  integer :: npts = 20, nj=20

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer



  Pstar = powerAmpScalar

  call delete_file('ai_predic.dat')
  call delete_file('ai_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp

  mumin=(10._kp)**(-3)
  mumax=0.512378_kp*0.99_kp
  
  do j=0,nj
    mu=mumin*(mumax/mumin)**(real(j,kp)/real(nj,kp))


  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = ai_lnrhoend(mu,Pstar)

  print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  print*,'mu=',mu,'xEnd=',ai_x_endinf(mu)

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

	xstar = ai_x_star(mu,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',ai_epsilon_one(xstar,mu)
 

       eps1 = ai_epsilon_one(xstar,mu)
       eps2 = ai_epsilon_two(xstar,mu)
       eps3 = ai_epsilon_three(xstar,mu)


       logErehGeV = log_energy_reheat_ingev(lnRhoReh)


       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('ai_predic.dat',mu,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('ai_nsr.dat',mu,ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do


end program aimain
