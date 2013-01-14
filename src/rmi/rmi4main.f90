!test the reheating derivation from slow-roll
program rmi4main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rmi4sr, only : rmi4_epsilon_one, rmi4_epsilon_two, rmi4_epsilon_three
  use rmi4reheat, only : rmi4_lnrhoend, rmi4_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,l
  integer :: npts = 20

  integer :: Nc, Nphi0, Nxend
  real(kp) ::cmin, cmax, phi0min, phi0max, xendmin, xendmax, c, phi0, xend

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!    Slow Roll Predictions   !!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Pstar = powerAmpScalar

  call delete_file('rmi4_predic.dat')
  call delete_file('rmi4_nsr.dat')

 Nc=10
 Nphi0=20
 Nxend=100

!  w = 1._kp/3._kp
  w=0._kp

 cmin=10._kp**(-3._kp)
 cmax=10._kp**(0._kp)

! do j=0,3 

! if (j .eq. 0) then
! c=-10._kp**(-8._kp) 
! end if
! if (j .eq. 1) then
! c=-10._kp**(-5._kp) 
! end if
! if (j .eq. 2) then
! c=-10._kp**(-4._kp) 
! end if
! if (j .eq. 3) then 
! c=-10._kp**(-3._kp) 
! end if

 !c=-10._kp**(-5._kp) 
 !c=-10._kp**(-4._kp) 
 c=-10._kp**(-2._kp) 

   phi0max=1._kp/sqrt(abs(c))
   phi0min=phi0max/(10._kp**4.)

    do k=0,Nphi0 
      phi0=phi0min*(phi0max/phi0min)**(real(k,kp)/Nphi0)  !logarithmic step

    xendmin =  1._kp+sqrt(2._kp*epsilon(1._kp)*(1._kp+c*phi0**2/4._kp)**2/(c**2*phi0**2)) !Using an asymptotic expression for eps1 when x->1, and requiring eps1>epsilon(1._kp) for numerical convergence
    xendmax = exp(1._kp)

    do l=0,Nxend 
      xend=xendmin*(xendmax/xendmin)**(real(l,kp)/Nxend)  !logarithmic step

      lnRhoRehMin = lnRhoNuc
      lnRhoRehMax = rmi4_lnrhoend(c,phi0,xend,Pstar)


      print *,'c=',c,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

      do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = rmi4_x_star(c,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar
 

       eps1 = rmi4_epsilon_one(xstar,c,phi0)
       eps2 = rmi4_epsilon_two(xstar,c,phi0)
       eps3 = rmi4_epsilon_three(xstar,c,phi0)
   

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)

       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('rmi4_predic.dat',c,phi0,xend,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('rmi4_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
      end do

    end do

  end do

! end do


end program rmi4main
