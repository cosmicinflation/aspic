!test the reheating derivation from slow-roll
program hytest
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use hyreheat, only : hy_x_reheat, hy_lnrhoend
  use hysrevol, only : hy_x_endinf, hy_epsilon_one, hy_epsilon_two
  use infinout, only : delete_file, livewrite
  implicit none
 
  real(kp) :: Pstar

  integer :: i
  integer :: npts = 20

  real(kp) :: p,mu,xstop,wreh,bfold
  real(kp) :: lnRhoReh,phi,eps1,eps2,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  p = 5._kp 
  mu = 0.6_kp  
  xstop = 1e-3

  wreh = 0.
 
  Pstar = powerAmpScalar

  call delete_file('hybrid_eps2eps1.dat')
  call delete_file('hybrid_nsr.dat')

 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = hy_lnrhoend(p,mu,xstop,Pstar)
  
  print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     phi = hy_x_reheat(p,mu,xstop,wreh,lnRhoReh,Pstar,bfold)

     print *,'lnRhoReh bfold= ',lnRhoReh,bfold

     eps1 = hy_epsilon_one(phi,p,mu)
     eps2 = hy_epsilon_two(phi,p,mu)

     call livewrite('hybrid_eps2eps1.dat',eps2,eps1,abs(bfold),lnRhoReh)

     ns = 1._kp - 2._kp*eps1 - eps2
     r =16._kp*eps1

     call livewrite('hybrid_nsr.dat',ns,r,abs(bfold),lnRhoReh)
  
  end do

  

end program hytest
