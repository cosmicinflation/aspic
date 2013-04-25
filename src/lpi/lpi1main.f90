!test the reheating derivation from slow-roll
program lpi1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lpi1sr, only : lpi1_epsilon_one, lpi1_epsilon_two, lpi1_epsilon_three
  use lpi1reheat, only : lpi1_lnrhoreh_max, lpi1_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev


  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nphi0

  real(kp) :: phi0,p,q,w,bfoldstar,phi0min,phi0max,pmin,pmax,qmin,qmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  integer :: npq
  real(kp), dimension(10) :: pvalues, qvalues

  real(kp) ::x,xmin,xmax


  Pstar = powerAmpScalar



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('lpi1_predic.dat')
  call delete_file('lpi1_nsr.dat')

  w=0._kp
!  w = 1._kp/3._kp

  npts = 20

  npq=3
 
  pvalues(1)=4._kp
  qvalues(1)=2._kp

  pvalues(2)=4._kp
  qvalues(2)=1._kp

  pvalues(3)=4._kp
  qvalues(3)=3._kp
  

  phi0min=10._kp**(-3._kp)
  phi0max=10._kp**(3._kp)

 do k=1,3

    p=pvalues(k)
    q=qvalues(k)

 if (k .eq. 1) nphi0=50
 if (k .eq. 2) nphi0=100
 if (k .eq. 3) nphi0=100


  do j=0,nphi0
    phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(nphi0,kp))
 
    lnRhoRehMin = lnRhoNuc
    lnRhoRehMax = lpi1_lnrhoreh_max(p,q,phi0,Pstar)


    print *,'phi0=',phi0,'p=',p,'q=',q,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

    do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = lpi1_x_star(p,q,phi0,w,lnRhoReh,Pstar,bfoldstar)


       eps1 = lpi1_epsilon_one(xstar,p,q,phi0)
       eps2 = lpi1_epsilon_two(xstar,p,q,phi0)
       eps3 = lpi1_epsilon_three(xstar,p,q,phi0)


       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('lpi1_predic.dat',p,q,phi0,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('lpi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

enddo


end program lpi1main
