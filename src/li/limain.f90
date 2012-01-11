!test the reheating derivation from slow-roll
program limain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lisr, only : li_epsilon_one, li_epsilon_two, li_epsilon_three
  use lireheat, only : li_lnrhoend, li_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 20

  real(kp) :: alpha,mu,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer
  
  real(kp) ::alphamin,alphamax,mumin,mumax

  Pstar = powerAmpScalar

  call delete_file('li_predic.dat')
  call delete_file('li_nsr.dat')

  alphamin=0.002
  alphamax=100000._kp
  mumin=0.0000001
  mumax=100._kp

!  w = 1._kp/3._kp
  w=0._kp

  do j=0,20
    mu=mumin*(mumax/mumin)**(real(j,kp)/real(20,kp))
    do k=0,20
      alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(20,kp))

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = li_lnrhoend(alpha,mu,Pstar)

  print *,'mu=',mu,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = li_x_star(alpha,mu,w,lnRhoReh,Pstar,bfoldstar)

       !print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

       eps1 = li_epsilon_one(xstar,alpha,mu)
       eps2 = li_epsilon_two(xstar,alpha,mu)
       eps3 = li_epsilon_three(xstar,alpha,mu)

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('li_predic.dat',alpha,mu,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('li_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

  end do

 end do



end program limain
