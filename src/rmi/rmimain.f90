!test the reheating derivation from slow-roll
program rmtest
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rmreheat, only : rm_x_star, rm_lnrhoend
  use rmsrevol, only : rm_x_endinf, rm_epsilon_one, rm_epsilon_two
  use infinout, only : delete_file, livewrite
  implicit none
 
  real(kp) :: Pstar

  integer :: i
  integer :: npts = 20

  real(kp) :: p,mu,nu,xstop,wreh,bfold
  real(kp) :: lnRhoReh,phi,eps1,eps2,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  p = 2._kp 
  mu = 0.5_kp
  nu = 0.01_kp
  xstop = 2._kp

  wreh = 0.
 
  Pstar = powerAmpScalar

  call delete_file('rmass_eps2eps1.dat')
  call delete_file('rmass_nsr.dat')

 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = rm_lnrhoend(p,mu,nu,xstop,Pstar)

  print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     phi = rm_x_star(p,mu,nu,xstop,wreh,lnRhoReh,Pstar,bfold)

     print *,'lnRhoReh bfold= ',lnRhoReh,bfold

     eps1 = rm_epsilon_one(phi,p,mu,nu)
     eps2 = rm_epsilon_two(phi,p,mu,nu)

     call livewrite('rmass_eps2eps1.dat',eps2,eps1,abs(bfold),lnRhoReh)

     ns = 1._kp - 2._kp*eps1 - eps2
     r =16._kp*eps1

     call livewrite('rmass_nsr.dat',ns,r,abs(bfold),lnRhoReh)
  
  end do

  

end program rmtest
