!test the reheating derivation from slow-roll
program ssbi4main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ssbi4sr, only : ssbi4_epsilon_one, ssbi4_epsilon_two, ssbi4_epsilon_three
  use ssbi4reheat, only : ssbi4_lnrhoend, ssbi4_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 20

  integer :: Nalpha,Nbeta
  real(kp) ::alphamin, alphamax, betamin, betamax, alpha, beta

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer
  real(kp), dimension(:), allocatable :: betavalues



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!    Slow Roll Predictions   !!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Pstar = powerAmpScalar

  call delete_file('ssbi4_predic.dat')
  call delete_file('ssbi4_nsr.dat')

 Nalpha=100
!  w = 1._kp/3._kp
  w=0._kp

  Nbeta=3
  allocate(betavalues(1:Nbeta))
  betavalues(1)=-10._kp**(-5._kp) 
  betavalues(2)=-10._kp**(-4._kp)
  betavalues(3)=-10._kp**(-3._kp)

    do j=1,Nbeta
    beta=betavalues(j)

   alphamin=10._kp**(-4._kp)
   alphamax=10._kp**(1._kp)

    do k=0,Nalpha 
      alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/Nalpha)  !logarithmic step

      lnRhoRehMin = lnRhoNuc
      lnRhoRehMax = ssbi4_lnrhoend(alpha,beta,Pstar)


      print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

      do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = ssbi4_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar
 

       eps1 = ssbi4_epsilon_one(xstar,alpha,beta)
       eps2 = ssbi4_epsilon_two(xstar,alpha,beta)
       eps3 = ssbi4_epsilon_three(xstar,alpha,beta)
   

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)

       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('ssbi4_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('ssbi4_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
      end do

    end do

 end do


end program ssbi4main
