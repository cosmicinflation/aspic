!test the reheating derivation from slow-roll
program mlfimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use mlfisr, only : mlfi_epsilon_one, mlfi_epsilon_two, mlfi_epsilon_three
  use mlfireheat, only : mlfi_lnrhoend 
  use mlfireheat, only : mlfi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  implicit none
 

  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i
  integer :: npts = 20

  real(kp) :: p,q,alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  Pstar = powerAmpScalar

!  call delete_file('mlfi_predic.dat')
!  call delete_file('mlfi_nsr.dat')

  p = 2._kp 
  q = 2._kp
  alpha =2._kp*(10._kp)**(-0._kp)
  w = 0._kp
 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = mlfi_lnrhoend(p,q,alpha,Pstar)

  print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = mlfi_x_star(p,q,alpha,w,lnRhoReh,Pstar,bfoldstar)

     print *,'lnRhoReh',lnRhoReh,'bfoldstar= ',bfoldstar

     eps1 = mlfi_epsilon_one(xstar,p,q,alpha)
     eps2 = mlfi_epsilon_two(xstar,p,q,alpha)
     eps3 = mlfi_epsilon_three(xstar,p,q,alpha)

     logErehGeV = log_energy_reheat_ingev(lnRhoReh)
     Treh =  10._kp**(logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp))


     ns = 1._kp - 2._kp*eps1 - eps2
     r =16._kp*eps1

     call livewrite('mlfi_predic.dat',alpha,p,q,eps1,eps2,eps3,r,ns,Treh)

     call livewrite('mlfi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
  end do

  

end program mlfimain
