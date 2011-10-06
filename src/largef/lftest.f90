!test the reheating derivation from slow-roll
program lftest
  use infprec, only : kp, transfert
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lfsrevol, only : lf_epsilon_one, lf_epsilon_two
  use lfreheat, only : lf_lnrhoend, lf_lnrhoreh 
  use lfreheat, only : lf_x_obs, lf_x_reheat
  use infinout, only : delete_file, livewrite
  implicit none

  type(transfert) :: lfData
  real(kp) :: Pstar, logErehGeV

  integer :: i
  integer :: npts = 20

  real(kp) :: p,wreh,bfold
  real(kp) :: lnRhoReh,phi,eps1,eps2,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  p = 2_kp 
  wreh = (p-2)/(p+2)

  lfData%real1 = p
  lfData%real2 = wreh


  Pstar = powerAmpScalar

  call delete_file('lfreh_eps2eps1.dat')
  call delete_file('lfreh_nsr.dat')

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = lf_lnrhoend(p,Pstar)

  print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     phi = lf_x_reheat(p,wreh,lnRhoReh,Pstar,bfold)

     print *,'lnRhoReh bfold= ',lnRhoReh,bfold

     eps1 = lf_epsilon_one(phi,p)
     eps2 = lf_epsilon_two(phi,p)

     logErehGeV = 0.25_kp*(lnRhoReh + 169.34_kp)/log(10._kp)
!     call livewrite('lfreh_eps2eps1.dat',eps2,eps1,abs(bfold),lnRhoReh)
     call livewrite('lfreh_eps2eps1.dat',eps2,eps1,abs(bfold),logErehGeV,lnRhoReh)

     ns = 1._kp - 2._kp*eps1 - eps2
     r =16._kp*eps1

     call livewrite('lfreh_nsr.dat',ns,r,abs(bfold),lnRhoReh)
  
  end do

  

end program lftest
