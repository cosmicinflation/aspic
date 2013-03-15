!test the reheating derivation from slow-roll
program sbimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use sbisr, only : sbi_epsilon_one, sbi_epsilon_two, sbi_epsilon_three, sbi_alpha_min
  use sbireheat, only : sbi_lnrhoend, sbi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev


  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nbeta

  real(kp) :: alpha,beta,w,bfoldstar,alphamin,alphamax,betamin,betamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer


  Pstar = powerAmpScalar



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    w=0._kp
!  w = 1._kp/3._kp

  call delete_file('sbi_predic.dat')
  call delete_file('sbi_nsr.dat')

  npts = 20

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    beta=0.00005    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  beta=0.00005
  alphamin=sbi_alpha_min(beta)*(1._kp+epsilon(1._kp))
  alphamax=1000._kp*alphamin
  alphamax=60._kp*alphamin
  nalpha=20

  do k=0,nalpha
       alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !log step

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = sbi_lnrhoend(alpha,beta,Pstar)

  print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = sbi_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)


       eps1 = sbi_epsilon_one(xstar,alpha,beta)
       eps2 = sbi_epsilon_two(xstar,alpha,beta)
       eps3 = sbi_epsilon_three(xstar,alpha,beta)


       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('sbi_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('sbi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!     beta=0.001     !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  beta=0.001
  alphamin=sbi_alpha_min(beta)*(1._kp+epsilon(1._kp))
  alphamax=20._kp*alphamin
  nalpha=30


  do k=0,nalpha
       alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !log step

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = sbi_lnrhoend(alpha,beta,Pstar)

  print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = sbi_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)


       eps1 = sbi_epsilon_one(xstar,alpha,beta)
       eps2 = sbi_epsilon_two(xstar,alpha,beta)
       eps3 = sbi_epsilon_three(xstar,alpha,beta)


       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('sbi_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('sbi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

end program sbimain
