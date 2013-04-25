!test the reheating derivation from slow-roll
program ssbi6main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ssbi6sr, only : ssbi6_epsilon_one, ssbi6_epsilon_two, ssbi6_epsilon_three, ssbi6_alphamax
  use ssbi6reheat, only : ssbi6_lnrhoreh_max, ssbi6_x_star
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

 call delete_file('ssbi6_abs_alpha_min.dat')

 Nbeta=10000
 betamin=0.000001_kp
 betamax=100._kp

 do j=0,Nbeta 
  beta=betamin*(betamax/betamin)**(real(j,kp)/Nbeta)  !logarithmic step
!  beta=betamin+(betamax-betamin)*(real(j,kp)/Nbeta)  !arithmetic step
  alpha=ssbi6_alphamax(beta)
  call livewrite('ssbi6_alphamax.dat',beta,alpha)
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

  call delete_file('ssbi6_predic.dat')
  call delete_file('ssbi6_nsr.dat')

 Nalpha=100

!  w = 1._kp/3._kp
  w=0._kp

  Nbeta=3
  allocate(betavalues(1:Nbeta))
  betavalues(1)=10._kp**(-5._kp)
  betavalues(2)=0.1_kp
  betavalues(3)=10._kp**(0._kp)

    do j=1,Nbeta
    beta=betavalues(j)

   alphamin=ssbi6_alphamax(beta)*1.001_kp
   if (alphamin .eq. 0._kp) alphamin=-10._kp**(-4._kp)
   alphamax=alphamin*10._kp**(4._kp)
   if (beta .eq. 10._kp**(-6._kp))   alphamax=alphamin*10._kp**(1.5_kp)

    do k=0,Nalpha 
      alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/Nalpha)  !logarithmic step

      lnRhoRehMin = lnRhoNuc
      lnRhoRehMax = ssbi6_lnrhoreh_max(alpha,beta,Pstar)


      print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

      do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = ssbi6_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar
 

       eps1 = ssbi6_epsilon_one(xstar,alpha,beta)
       eps2 = ssbi6_epsilon_two(xstar,alpha,beta)
       eps3 = ssbi6_epsilon_three(xstar,alpha,beta)
   

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)

       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('ssbi6_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('ssbi6_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
      end do

    end do

 end do


end program ssbi6main
