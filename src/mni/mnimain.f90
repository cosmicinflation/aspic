!test the reheating derivation from slow-roll
program mnimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use mnisr, only : mni_epsilon_one, mni_epsilon_two, mni_epsilon_three
  use mnireheat, only : mni_lnrhoend 
  use mnireheat, only : mni_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  implicit none
 

  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i
  integer :: npts = 20

  real(kp) :: f,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  Pstar = powerAmpScalar

!  call delete_file('mni_predic.dat')
!  call delete_file('mni_nsr.dat')

  f =100._kp*10._kp**(0)
  w = 0._kp

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = mni_lnrhoend(f,Pstar)

  print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = mni_x_star(f,w,lnRhoReh,Pstar,bfoldstar)

     print *,'lnRhoReh',lnRhoReh,'bfoldstar= ',bfoldstar

     eps1 = mni_epsilon_one(xstar,f)
     eps2 = mni_epsilon_two(xstar,f)
     eps3 = mni_epsilon_three(xstar,f)

     logErehGeV = log_energy_reheat_ingev(lnRhoReh)
     Treh =  10._kp**(logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp))


     ns = 1._kp - 2._kp*eps1 - eps2
     r =16._kp*eps1

     call livewrite('mni_predic.dat',f,eps1,eps2,eps3,r,ns,Treh)

     call livewrite('mni_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
  end do

  

end program mnimain
