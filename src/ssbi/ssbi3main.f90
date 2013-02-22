!test the reheating derivation from slow-roll
program ssbi3main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ssbi3sr, only : ssbi3_epsilon_one, ssbi3_epsilon_two, ssbi3_epsilon_three, ssbi3_alphamin
  use ssbi3reheat, only : ssbi3_lnrhoend, ssbi3_x_star
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
!!!!!!!         Prior Space        !!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call delete_file('ssbi3_alphamin.dat')

 Nbeta=1000
 betamin=0.000001_kp!10._kp**(-3._kp)
 betamax=100._kp!10._kp**(3._kp)

 do j=0,Nbeta 
  beta=-betamin*(betamax/betamin)**(real(j,kp)/Nbeta)  !logarithmic step
!  beta=betamin+(betamax-betamin)*(real(j,kp)/Nbeta)  !arithmetic step
  alpha=ssbi3_alphamin(beta)
  call livewrite('ssbi3_alphamin.dat',beta,alpha)
 end do

print *,'Priors Written'




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!    Slow Roll Predictions   !!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Pstar = powerAmpScalar

  call delete_file('ssbi3_predic.dat')
  call delete_file('ssbi3_nsr.dat')

 Nalpha=100
 Nbeta=3

!  w = 1._kp/3._kp
  w=0._kp

  Nbeta=3
  allocate(betavalues(1:Nbeta))
  betavalues(1)=-10._kp**(-3._kp)
  betavalues(2)=-5.*10._kp**(-3._kp)
  betavalues(3)=-10._kp**(-2._kp)

    do j=1,Nbeta
    beta=betavalues(j)

   alphamin=ssbi3_alphamin(beta)*(1._kp+1000._kp*epsilon(1._kp))
   alphamax=alphamin*10._kp**(3._kp)

    do k=0,Nalpha 
      alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/Nalpha)  !logarithmic step

      lnRhoRehMin = lnRhoNuc
      lnRhoRehMax = ssbi3_lnrhoend(alpha,beta,Pstar)


      print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

      do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = ssbi3_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar
 

       eps1 = ssbi3_epsilon_one(xstar,alpha,beta)
       eps2 = ssbi3_epsilon_two(xstar,alpha,beta)
       eps3 = ssbi3_epsilon_three(xstar,alpha,beta)
   

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)

       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('ssbi3_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('ssbi3_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
      end do

    end do

 end do


end program ssbi3main
