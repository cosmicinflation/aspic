!test the reheating derivation from slow-roll
program lmi2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lmi2sr, only : lmi2_epsilon_one, lmi2_epsilon_two, lmi2_epsilon_three, lmi2_xin_min
  use lmi2reheat, only : lmi2_lnrhoend, lmi2_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
!  use cosmopar, only : QrmsOverT

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,NxEnd
  integer :: npts = 20
         
  real(kp), dimension(1:6) :: gamma_lmi_values

  integer, dimension(1:6) :: NxEnd_values              

  real(kp) :: xEndMin              !to be specified by lmi2_xin_min
  real(kp) :: xEndMax    

  integer :: Nbeta=10
  real(kp) :: betamin=0.1
  real(kp) :: betamax=10.

  real(kp) :: gamma_lmi,beta,xEnd,w,bfoldstar,alpha
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
gamma_lmi_values(1)=0.45_kp
NxEnd_values(1)=100
gamma_lmi_values(2)=0.5_kp
NxEnd_values(2)=50
gamma_lmi_values(3)=0.55_kp
NxEnd_values(3)=30
gamma_lmi_values(4)=0.6_kp
NxEnd_values(4)=20
gamma_lmi_values(5)=0.63_kp
NxEnd_values(5)=15
gamma_lmi_values(6)=0.68_kp
NxEnd_values(6)=15
endif

if (beta.eq.10) then
gamma_lmi_values(1)=0.18_kp
NxEnd_values(1)=100
gamma_lmi_values(2)=0.2_kp
NxEnd_values(2)=50
gamma_lmi_values(3)=0.22_kp
NxEnd_values(3)=20
gamma_lmi_values(4)=0.235_kp
NxEnd_values(4)=20
gamma_lmi_values(5)=0.25_kp
NxEnd_values(5)=20
gamma_lmi_values(6)=0.27_kp
NxEnd_values(6)=20
endif

if (beta.eq.0.01) then
gamma_lmi_values(1)=0.1_kp
NxEnd_values(1)=300
gamma_lmi_values(2)=0.7_kp
NxEnd_values(2)=25
gamma_lmi_values(3)=0.85_kp
NxEnd_values(3)=20
gamma_lmi_values(4)=0.95_kp
NxEnd_values(4)=15
gamma_lmi_values(5)=0.985_kp
NxEnd_values(5)=20
gamma_lmi_values(6)=0.995_kp
NxEnd_values(6)=20
endif




 do j=1,size(gamma_lmi_values) 
 gamma_lmi=gamma_lmi_values(j)
 NxEnd=nxEnd_values(j)

alpha=4._kp*(1._kp-gamma_lmi)
xEndMin=lmi2_xin_min(gamma_lmi,beta)*1.1
xEndMax=100._kp*max(alpha,(beta*gamma_lmi)**(1._kp/(1._kp-gamma_lmi)),(alpha*beta*gamma_lmi)**(1._kp/(2._kp-gamma_lmi)))


if (beta.eq.0.01) then
xEndMax=1._kp*max(sqrt(alpha),alpha,(beta*gamma_lmi)**(1._kp/(1._kp-gamma_lmi)), &
        (alpha*beta*gamma_lmi)**(1._kp/(2._kp-gamma_lmi)))
endif

 do k=0,NxEnd
 xEnd=xEndMin*(xEndMax/xEndMin)**(real(k,kp)/NxEnd)  !logarithmic step
! xEnd=xEndMin+(xEndMax-xEndMin)*(real(k,kp)/NxEnd)  !arithmetic step

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = lmi2_lnrhoend(gamma_lmi,beta,xEnd,Pstar)

  print *,'gamma_lmi=',gamma_lmi,'beta=',beta,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)


       xstar = lmi2_x_star(gamma_lmi,beta,xEnd,w,lnRhoReh,Pstar,bfoldstar)



       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar


       eps1 = lmi2_epsilon_one(xstar,gamma_lmi,beta)
       eps2 = lmi2_epsilon_two(xstar,gamma_lmi,beta)
       eps3 = lmi2_epsilon_three(xstar,gamma_lmi,beta)
   

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)

       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('lmi2_predic.dat',gamma_lmi,beta,xEnd,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('lmi2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

 end do


 

end program lmi2main
