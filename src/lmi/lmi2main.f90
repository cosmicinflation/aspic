!test the reheating derivation from slow-roll
program lmi2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lmi2sr, only : lmi2_epsilon_one, lmi2_epsilon_two, lmi2_epsilon_three, lmi2_xini_min
  use lmi2reheat, only : lmi2_lnrhoend, lmi2_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
!  use cosmopar, only : QrmsOverT

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,NxEnd
  integer :: npts = 20
         
  real(kp), dimension(1:6) :: gamValues

  integer, dimension(1:6) :: NxEndValues              

  real(kp) :: xEndMin              !to be specified by lmi2_xini_min
  real(kp) :: xEndMax    

  integer :: Nbeta=10
  real(kp) :: betamin=0.1
  real(kp) :: betamax=10.

  real(kp) :: gam,beta,xEnd,w,bfoldstar,alpha
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer



  Pstar = powerAmpScalar

  call delete_file('lmi2_predic.dat')
  call delete_file('lmi2_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp

beta=0.01
beta=1.
!beta=10.


if (beta.eq.1) then
gamValues(1)=0.45_kp
NxEndValues(1)=100
gamValues(2)=0.5_kp
NxEndValues(2)=50
gamValues(3)=0.55_kp
NxEndValues(3)=30
gamValues(4)=0.6_kp
NxEndValues(4)=20
gamValues(5)=0.63_kp
NxEndValues(5)=15
gamValues(6)=0.68_kp
NxEndValues(6)=15
endif

if (beta.eq.10) then
gamValues(1)=0.18_kp
NxEndValues(1)=100
gamValues(2)=0.2_kp
NxEndValues(2)=50
gamValues(3)=0.22_kp
NxEndValues(3)=20
gamValues(4)=0.235_kp
NxEndValues(4)=20
gamValues(5)=0.25_kp
NxEndValues(5)=20
gamValues(6)=0.27_kp
NxEndValues(6)=20
endif

if (beta.eq.0.01) then
gamValues(1)=0.1_kp
NxEndValues(1)=300
gamValues(2)=0.7_kp
NxEndValues(2)=25
gamValues(3)=0.85_kp
NxEndValues(3)=20
gamValues(4)=0.95_kp
NxEndValues(4)=15
gamValues(5)=0.985_kp
NxEndValues(5)=20
gamValues(6)=0.995_kp
NxEndValues(6)=20
endif




 do j=1,size(gamValues) 
 gam=gamValues(j)
 NxEnd=nxEndValues(j)

alpha=4._kp*(1._kp-gam)
xEndMin=lmi2_xini_min(gam,beta)*1.1
xEndMax=100._kp*max(alpha,(beta*gam)**(1._kp/(1._kp-gam)),(alpha*beta*gam)**(1._kp/(2._kp-gam)))


if (beta.eq.0.01) then
xEndMax=1._kp*max(sqrt(alpha),alpha,(beta*gam)**(1._kp/(1._kp-gam)), &
        (alpha*beta*gam)**(1._kp/(2._kp-gam)))
endif

 do k=0,NxEnd
 xEnd=xEndMin*(xEndMax/xEndMin)**(real(k,kp)/NxEnd)  !logarithmic step
! xEnd=xEndMin+(xEndMax-xEndMin)*(real(k,kp)/NxEnd)  !arithmetic step

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = lmi2_lnrhoend(gam,beta,xEnd,Pstar)

  print *,'gam=',gam,'beta=',beta,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)


       xstar = lmi2_x_star(gam,beta,xEnd,w,lnRhoReh,Pstar,bfoldstar)



       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar


       eps1 = lmi2_epsilon_one(xstar,gam,beta)
       eps2 = lmi2_epsilon_two(xstar,gam,beta)
       eps3 = lmi2_epsilon_three(xstar,gam,beta)
   

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)

       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('lmi2_predic.dat',gam,beta,xEnd,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('lmi2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

 end do


 

end program lmi2main
