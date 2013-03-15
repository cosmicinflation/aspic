!test the reheating derivation from slow-roll
program twimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use twisr, only : twi_epsilon_one, twi_epsilon_two, twi_epsilon_three
  use twireheat, only : twi_lnrhoend, twi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 20

  integer :: Nphi0=15
  real(kp) :: phi0min=0.001
  real(kp) :: phi0max=5.

  integer :: NxEnd=80
  real(kp) :: yEndmin=2.
  real(kp) :: yEndmax=20.

  real(kp) :: phi0,xEnd,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer



  Pstar = powerAmpScalar

  call delete_file('twi_predic.dat')
  call delete_file('twi_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp

 do j=0,Nphi0 
 phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/Nphi0)

 do k=0,NxEnd
 xEnd=phi0*yEndmin*(yEndmax/yEndmin)**(real(k,kp)/NxEnd)*(1.+epsilon(1._kp)) !logarithmic step
 xEnd=phi0*yEndmin+(yEndmax-yEndmin)*(real(k,kp)/NxEnd)*phi0*(1.+epsilon(1._kp)) !arithmetic step
 xEnd=phi0*yEndmin*(yEndmax/yEndmin)**(2.**(2.**(real(k,kp)/NxEnd)-1.)-1.)*(1.+epsilon(1._kp)) !ultralogarithmic step


  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = twi_lnrhoend(phi0,xEnd,Pstar)

  print *,'phi0=',phi0,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

	xstar = twi_x_star(phi0,xEnd,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',twi_epsilon_one(xstar,phi0)
 

       eps1 = twi_epsilon_one(xstar,phi0)
       eps2 = twi_epsilon_two(xstar,phi0)
       eps3 = twi_epsilon_three(xstar,phi0)


       logErehGeV = log_energy_reheat_ingev(lnRhoReh)


       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('twi_predic.dat',phi0,xEnd,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('twi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

end do


end program twimain
