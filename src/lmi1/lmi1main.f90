!test the reheating derivation from slow-roll
program lmi1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lmi1sr, only : lmi1_epsilon_one, lmi1_epsilon_two, lmi1_epsilon_three
  use lmi1reheat, only : lmi1_lnrhoend, lmi1_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
!  use cosmopar, only : QrmsOverT

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 20

  integer :: Ngamma_lmi=1000         !for beta=0.001: Ngamma=20
                                   !for beta=1: Ngamma=50
                                   !for beta=50: Ngamma=1000
  real(kp) :: gamma_lmimin=0.00005   !for beta = 0.001:  gammamin=0.004
                                   !for beta=1: gammamin=0.001
                                   !for beta=50: gammamin=0.00005
  real(kp) :: gamma_lmimax=0.1    !for beta=0.001: gammamax=0.99
                                   !for beta=1: gammamax=0.99
                                   !for beta=50: gammamax=0.1

  integer :: Nbeta=10
  real(kp) :: betamin=0.1
  real(kp) :: betamax=10.

  real(kp) :: gamma_lmi,alpha,beta,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer



  Pstar = powerAmpScalar

  call delete_file('lmi1_predic.dat')
  call delete_file('lmi1_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp

beta=0.001
beta=1.
beta=50.

 do j=0,Ngamma_lmi 
 gamma_lmi=gamma_lmimin*(gamma_lmimax/gamma_lmimin)**(real(j,kp)/Ngamma_lmi)  !logarithmic step
 gamma_lmi=gamma_lmimin+(gamma_lmimax-gamma_lmimin)*(real(j,kp)/Ngamma_lmi)  !arithmetic step
 gamma_lmi=sqrt(gamma_lmimin+(gamma_lmimax-gamma_lmimin)*(real(j,kp)/Ngamma_lmi))  !square root step



!  alpha=4.*(1.-gamma_lmi)
!  w=(alpha-2.)/(alpha+2.)

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = lmi1_lnrhoend(gamma_lmi,beta,Pstar)

  print *,'gamma_lmi=',gamma_lmi,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = lmi1_x_star(gamma_lmi,beta,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar
 

       eps1 = lmi1_epsilon_one(xstar,gamma_lmi,beta)
       eps2 = lmi1_epsilon_two(xstar,gamma_lmi,beta)
       eps3 = lmi1_epsilon_three(xstar,gamma_lmi,beta)
   

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)

       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('lmi1_predic.dat',gamma_lmi,beta,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('lmi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do


 

end program lmi1main
