!test the reheating derivation from slow-roll
program ostimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ostisr, only : osti_epsilon_one, osti_epsilon_two, osti_epsilon_three
  use ostireheat, only : osti_lnrhoreh_max, osti_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  integer :: Nphi0
  real(kp) :: phi0min=4._kp
  real(kp) :: phi0max=10._kp**3

  real(kp) :: phi0,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  Nphi0=40

  Pstar = powerAmpScalar

  call delete_file('osti_predic.dat')
  call delete_file('osti_nsr.dat')

!  w = 1._kp/3._kp
  w=0._kp

 do j=0,Nphi0+2

 phi0=phi0min+(phi0max-phi0min)*(real(j,kp)/real(Nphi0,kp)) !arithmetic step
 phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(Nphi0,kp)) !logarithmic step

   if (j .eq. Nphi0+1) phi0=10._kp**4.
   if (j .eq. Nphi0+2) phi0=10._kp**5.

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = osti_lnrhoreh_max(phi0,Pstar)

  print *,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

	xstar = osti_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar
 
       eps1 = osti_epsilon_one(xstar,phi0)
       eps2 = osti_epsilon_two(xstar,phi0)
       eps3 = osti_epsilon_three(xstar,phi0)


       logErehGeV = log_energy_reheat_ingev(lnRhoReh)


       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('osti_predic.dat',phi0,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('osti_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

end program ostimain
