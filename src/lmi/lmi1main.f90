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

  integer :: Ngam                  !for beta=0.001: Ngam=20
                                   !for beta=1: Ngam=50
                                   !for beta=50: Ngam=1000
  real(kp) :: gammin               !for beta = 0.001:  gammin=0.004
                                   !for beta=1: gammin=0.001
                                   !for beta=50: gammin=0.00005
  real(kp) :: gammax               !for beta=0.001: gammax=0.99
                                   !for beta=1: gammax=0.99
                                   !for beta=50: gammax=0.1

  integer :: Nbeta=10
  real(kp) :: betamin=0.1
  real(kp) :: betamax=10.

  real(kp) :: gam,alpha,beta,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer



  Pstar = powerAmpScalar

  call delete_file('lmi1_predic.dat')
  call delete_file('lmi1_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp

!!!!!!!!!!!!!!!!!!!!!!!
!!!!!  beta=0.001 !!!!!
!!!!!!!!!!!!!!!!!!!!!!!

  beta=0.001
  Ngam=20
  gammin=0.004
  gammax=0.99


 do j=0,Ngam 
 gam=gammin*(gammax/gammin)**(real(j,kp)/Ngam)  !logarithmic step
 gam=gammin+(gammax-gammin)*(real(j,kp)/Ngam)  !arithmetic step
 gam=sqrt(gammin+(gammax-gammin)*(real(j,kp)/Ngam))  !square root step

!  alpha=4.*(1.-gam)
!  w=(alpha-2.)/(alpha+2.)

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = lmi1_lnrhoend(gam,beta,Pstar)

  print *,'gam=',gam,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = lmi1_x_star(gam,beta,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar
 

       eps1 = lmi1_epsilon_one(xstar,gam,beta)
       eps2 = lmi1_epsilon_two(xstar,gam,beta)
       eps3 = lmi1_epsilon_three(xstar,gam,beta)
   

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)

       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('lmi1_predic.dat',gam,beta,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('lmi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

!!!!!!!!!!!!!!!!!!!!!!!
!!!!!    beta=1   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!

  beta=1.
  gammin=0.001
  gammax=0.62
  Ngam=30

 do j=0,Ngam 
 gam=gammin*(gammax/gammin)**(real(j,kp)/Ngam)  !logarithmic step
 gam=gammin+(gammax-gammin)*(real(j,kp)/Ngam)  !arithmetic step
 gam=sqrt(gammin+(gammax-gammin)*(real(j,kp)/Ngam))  !square root step


!  alpha=4.*(1.-gam)
!  w=(alpha-2.)/(alpha+2.)

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = lmi1_lnrhoend(gam,beta,Pstar)

  print *,'gam=',gam,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = lmi1_x_star(gam,beta,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar
 

       eps1 = lmi1_epsilon_one(xstar,gam,beta)
       eps2 = lmi1_epsilon_two(xstar,gam,beta)
       eps3 = lmi1_epsilon_three(xstar,gam,beta)
   

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)

       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('lmi1_predic.dat',gam,beta,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('lmi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do


!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   beta=50   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!
  beta=50.
  gammin=0.005
  gammax=0.07
  Ngam=50

 do j=0,Ngam 
 gam=gammin*(gammax/gammin)**(real(j,kp)/Ngam)  !logarithmic step
 gam=gammin+(gammax-gammin)*(real(j,kp)/Ngam)  !arithmetic step
 gam=gammin+(gammax-gammin)*sqrt(real(j,kp)/real(Ngam,kp))  !square root step


!  alpha=4.*(1.-gam)
!  w=(alpha-2.)/(alpha+2.)

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = lmi1_lnrhoend(gam,beta,Pstar)

  print *,'gam=',gam,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = lmi1_x_star(gam,beta,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar
 

       eps1 = lmi1_epsilon_one(xstar,gam,beta)
       eps2 = lmi1_epsilon_two(xstar,gam,beta)
       eps3 = lmi1_epsilon_three(xstar,gam,beta)
   

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)

       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('lmi1_predic.dat',gam,beta,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('lmi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do



 

end program lmi1main
