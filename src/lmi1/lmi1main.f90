!test the reheating derivation from slow-roll
program lmi1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lmi1sr, only : lmi1_epsilon_one, lmi1_epsilon_two, lmi1_epsilon_three
  use lmi1sr, only : lmi1_beta,lmi1_alpha
  use lmi1reheat, only : lmi1_lnrhoend, lmi1_x_star,lmi1_M
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
!  use cosmopar, only : QrmsOverT

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  integer :: Ngamma_lmi=10
  real(kp) :: gamma_lmimin=0.7
  real(kp) :: gamma_lmimax=0.99

  real(kp) :: gamma_lmi,M,alpha,beta,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer



  Pstar = powerAmpScalar

  call delete_file('lmi1_predic.dat')
  call delete_file('lmi1_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp

 do j=0,Ngamma_lmi 
 gamma_lmi=gamma_lmimin*(gamma_lmimax/gamma_lmimin)**(real(j,kp)/Ngamma_lmi)
 M=20._kp

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = lmi1_lnrhoend(gamma_lmi,M,Pstar)

  print *,'gamma_lmi=',gamma_lmi,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       M=lmi1_M(gamma_lmi,lnRhoReh,w,Pstar)

       xstar = lmi1_x_star(gamma_lmi,M,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'M=',M,'xstar=',xstar
 

       eps1 = lmi1_epsilon_one(xstar,gamma_lmi,M)
       eps2 = lmi1_epsilon_two(xstar,gamma_lmi,M)
       eps3 = lmi1_epsilon_three(xstar,gamma_lmi,M)
   

       alpha=lmi1_alpha(gamma_lmi,M)
       beta=lmi1_beta(gamma_lmi,M)

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)


       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('lmi1_predic.dat',gamma_lmi,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('lmi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

 

end program lmi1main
