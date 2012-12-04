!test the reheating derivation from slow-roll
program ssi5main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ssi5sr, only : ssi5_epsilon_one, ssi5_epsilon_two, ssi5_epsilon_three, ssi5_abs_alpha_min
  use ssi5reheat, only : ssi5_lnrhoend, ssi5_x_star
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!         Prior Space        !!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 call delete_file('ssi5_alphamin.dat')

 Nbeta=1000
 betamin=0.0001_kp
 betamax=100._kp

 do j=0,Nbeta 
  beta=betamin*(betamax/betamin)**(real(j,kp)/Nbeta)  !logarithmic step
!  beta=betamin+(betamax-betamin)*(real(j,kp)/Nbeta)  !arithmetic step
  alpha=-ssi5_abs_alpha_min(beta)
  call livewrite('ssi5_abs_alpha_min.dat',beta,-alpha)
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

  call delete_file('ssi5_predic.dat')
  call delete_file('ssi5_nsr.dat')

 Nalpha=1000
 Nbeta=100

!  w = 1._kp/3._kp
  w=0._kp

 betamin=10._kp**(-3._kp)
 betamax=10._kp**(3._kp)

! do j=0,Nbeta 
! beta=betamin*(betamax/betamin)**(real(j,kp)/Nbeta)  !logarithmic step

  beta=10._kp**(-6._kp)
!  beta=10._kp**(-5._kp)
!  beta=10._kp**(-4._kp)

   alphamin=-ssi5_abs_alpha_min(beta)*1.001_kp
   alphamax=alphamin*10._kp**(1._kp)
   if (beta .eq. 10._kp**(-6._kp))   alphamax=alphamin*10._kp**(1.5_kp)

    do k=0,Nalpha 
      alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/Nalpha)  !logarithmic step

      lnRhoRehMin = lnRhoNuc
      lnRhoRehMax = ssi5_lnrhoend(alpha,beta,Pstar)


      print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

      do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = ssi5_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar
 

       eps1 = ssi5_epsilon_one(xstar,alpha,beta)
       eps2 = ssi5_epsilon_two(xstar,alpha,beta)
       eps3 = ssi5_epsilon_three(xstar,alpha,beta)
   

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)

       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('ssi5_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('ssi5_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
      end do

    end do

! end do


end program ssi5main
