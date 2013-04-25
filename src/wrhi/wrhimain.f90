!test the reheating derivation from slow-roll
program wrhimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use wrhisr, only : wrhi_epsilon_one, wrhi_epsilon_two, wrhi_epsilon_three
  use wrhireheat, only : wrhi_lnrhoreh_max, wrhi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

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

end program wrhimain
