!test the reheating derivation from slow-roll
program lmi3main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lmi3sr, only : lmi3_epsilon_one, lmi3_epsilon_two, lmi3_epsilon_three
  use lmi3sr, only : lmi3_beta,lmi3_alpha
  use lmi3reheat, only : lmi3_lnrhoend, lmi3_x_star,lmi3_M
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 20

  integer :: Ngamma_lmi=5
  real(kp) :: gamma_lmimin=0.1
  real(kp) :: gamma_lmimax=0.7

  integer :: NxEnd=5
  real(kp) :: xEndmin=100._kp
  real(kp) :: xEndmax=1000000._kp

  real(kp) :: gamma_lmi,M,xEnd,alpha,beta,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer



  Pstar = powerAmpScalar

  call delete_file('lmi3_predic.dat')
  call delete_file('lmi3_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp

 do j=0,Ngamma_lmi 
 gamma_lmi=gamma_lmimin*(gamma_lmimax/gamma_lmimin)**(real(j,kp)/Ngamma_lmi)

 do k=0,NxEnd 
 xEnd =xEndmin*(xEndmax/xEndmin)**(real(k,kp)/NxEnd )


 M=1._kp
! M=lmi3_M(gamma_lmi,xEnd,lnRhoNuc,w,Pstar)

 print*,'M=',M

 !pause

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = lmi3_lnrhoend(gamma_lmi,M,xEnd,Pstar)

  print *,'gamma_lmi=',gamma_lmi,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  !pause 

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       !M=lmi3_M(gamma_lmi,xEnd,lnRhoReh,w,Pstar)

       xstar = lmi3_x_star(gamma_lmi,M,xEnd,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'M=',M,'xstar=',xstar
 

       eps1 = lmi3_epsilon_one(xstar,gamma_lmi,M)
       eps2 = lmi3_epsilon_two(xstar,gamma_lmi,M)
       eps3 = lmi3_epsilon_three(xstar,gamma_lmi,M)
   

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)


       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('lmi3_predic.dat',gamma_lmi,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('lmi3_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

end do

 

end program lmi3main
