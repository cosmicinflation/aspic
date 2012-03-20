!test the reheating derivation from slow-roll
program lmi2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lmi2sr, only : lmi2_epsilon_one, lmi2_epsilon_two, lmi2_epsilon_three
  use lmi2sr, only : lmi2_beta,lmi2_alpha
  use lmi2reheat, only : lmi2_lnrhoend, lmi2_x_star,lmi2_M
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
!  use cosmopar, only : QrmsOverT

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  integer :: Ngamma_lmi=10
  real(kp) :: gamma_lmimin=0.95
  real(kp) :: gamma_lmimax=0.999

  real(kp) :: gamma_lmi,M,alpha,beta,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer



  Pstar = powerAmpScalar

  call delete_file('lmi2_predic.dat')
  call delete_file('lmi2_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp

 do j=0,Ngamma_lmi 
 gamma_lmi=gamma_lmimin*(gamma_lmimax/gamma_lmimin)**(real(j,kp)/Ngamma_lmi)
 M=0.0000006_kp
 !M=lmi2_M(gamma_lmi,lnRhoNuc,w,Pstar)
 print*,'M=',M

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = lmi2_lnrhoend(gamma_lmi,M,Pstar)

  print *,'gamma_lmi=',gamma_lmi,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  !pause 

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       !M=lmi2_M(gamma_lmi,lnRhoReh,w,Pstar)

       xstar = lmi2_x_star(gamma_lmi,M,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'M=',M,'xstar=',xstar
 

       eps1 = lmi2_epsilon_one(xstar,gamma_lmi,M)
       eps2 = lmi2_epsilon_two(xstar,gamma_lmi,M)
       eps3 = lmi2_epsilon_three(xstar,gamma_lmi,M)
   

       alpha=lmi2_alpha(gamma_lmi,M)
       beta=lmi2_beta(gamma_lmi,M)
       !print*, 'M=',720._kp*acos(-1._kp)**2*(alpha-beta*gamma_lmi*xstar**gamma_lmi)**2* &
       !         exp(beta*xstar**gamma_lmi)*xstar**(-alpha-2._kp)*QrmsOverT**2

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)


       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('lmi2_predic.dat',gamma_lmi,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('lmi2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

 

end program lmi2main
