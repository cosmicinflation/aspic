!test the reheating derivation from slow-roll
program mixlftest
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use mixlfreheat, only : mixlf_x_reheat, mixlf_lnrhoend
  use mixlfsrevol, only : mixlf_x_endinf, mixlf_epsilon_one, mixlf_epsilon_two
  use infinout, only : delete_file, livewrite
  implicit none
 
  real(kp) :: Pstar

  integer :: i
  integer :: npts = 20

  real(kp) :: p,q,alpha,wreh,bfold
  real(kp) :: lnRhoReh,phi,eps1,eps2,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  p = 2._kp 
  q = 4._kp
  alpha = 0.5_kp

  wreh = 0.
 
  Pstar = powerAmpScalar

  call delete_file('mixreh_eps2eps1.dat')
  call delete_file('mixreh_nsr.dat')

 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = mixlf_lnrhoend(p,q,alpha,Pstar)

  print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     phi = mixlf_x_reheat(p,q,alpha,wreh,lnRhoReh,Pstar,bfold)

     print *,'lnRhoReh bfold= ',lnRhoReh,bfold

     eps1 = mixlf_epsilon_one(phi,p,q,alpha)
     eps2 = mixlf_epsilon_two(phi,p,q,alpha)

     call livewrite('mixreh_eps2eps1.dat',eps2,eps1,abs(bfold),lnRhoReh)

     ns = 1._kp - 2._kp*eps1 - eps2
     r =16._kp*eps1

     call livewrite('mixreh_nsr.dat',ns,r,abs(bfold),lnRhoReh)
  
  end do

  

end program mixlftest
