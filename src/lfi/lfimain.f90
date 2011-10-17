!test the reheating derivation from slow-roll
program lfimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lfisr, only : lfi_epsilon_one, lfi_epsilon_two, lfi_epsilon_three
  use lfireheat, only : lfi_lnrhoend, lfi_lnrhoreh_fromepsilon 
  use lfireheat, only : lfi_xp_fromepsilon, lfi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  implicit none

  
  real(kp) :: Pstar, logErehGeV

  integer :: i
  integer :: npts = 20

  real(kp) :: p,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  p = 2_kp 
  w = (p-2)/(p+2)
 
  Pstar = powerAmpScalar

  call delete_file('lfi_predic.dat')
  call delete_file('lfi_nsr.dat')

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = lfi_lnrhoend(p,Pstar)

  print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = lfi_x_star(p,w,lnRhoReh,Pstar,bfoldstar)

     print *,'lnRhoReh bfoldstar= ',lnRhoReh,bfoldstar

     eps1 = lfi_epsilon_one(xstar,p)
     eps2 = lfi_epsilon_two(xstar,p)
     eps3 = lfi_epsilon_three(xstar,p)

     logErehGeV = log_energy_reheat_ingev(lnRhoReh)


     ns = 1._kp - 2._kp*eps1 - eps2
     r =16._kp*eps1

     call livewrite('lfi_predic.dat',p,eps1,eps2,eps3,r,ns,logErehGeV)

     call livewrite('lfi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
  end do

  

end program lfimain
