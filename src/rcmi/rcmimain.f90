!test the reheating derivation from slow-roll
program rcmimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rcmisr, only : rcmi_epsilon_one, rcmi_epsilon_two, rcmi_epsilon_three
  use rcmireheat, only : rcmi_lnrhoend 
  use rcmireheat, only : rcmi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  implicit none
 

  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(1:7) ::alphavalues

  alphavalues(1)=(10._kp)**(-6.)
  alphavalues(2)=5.*(10._kp)**(-5.)
  alphavalues(3)=(10._kp)**(-4.)
  alphavalues(4)=2.*(10._kp)**(-4.)
  alphavalues(5)=3.*(10._kp)**(-4.)
  alphavalues(6)=4.*(10._kp)**(-4.)
  alphavalues(7)=5.*(10._kp)**(-4.)

  Pstar = powerAmpScalar
  w = 0._kp

  call delete_file('rcmi_predic.dat')
  call delete_file('rcmi_nsr.dat')

  do j=1,size(alphavalues)
   
  alpha=alphavalues(j)

 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = rcmi_lnrhoend(alpha,Pstar)

  print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rcmi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

     print *,'lnRhoReh',lnRhoReh,'bfoldstar= ',bfoldstar

     eps1 = rcmi_epsilon_one(xstar,alpha)
     eps2 = rcmi_epsilon_two(xstar,alpha)
     eps3 = rcmi_epsilon_three(xstar,alpha)

     logErehGeV = log_energy_reheat_ingev(lnRhoReh)
     Treh =  10._kp**(logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp))

     ns = 1._kp - 2._kp*eps1 - eps2
     r =16._kp*eps1

     call livewrite('rcmi_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

     call livewrite('rcmi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
  end do

end do

  

end program rcmimain
