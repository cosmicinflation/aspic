!test the reheating derivation from slow-roll
program kksftest  
  use infprec, only : kp, transfert
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use kksfreheat, only : kksf_x_reheat, kksf_lnrhoend
  use kksfsrevol, only : kksf_x_endinf, kksf_epsilon_one, kksf_epsilon_two
  use infinout, only : delete_file, livewrite
  implicit none

  type(transfert) :: kksfData
  real(kp) :: Pstar,calF

  integer :: i
  integer :: npts = 20

  real(kp) :: mu,p,wreh,bfold
  real(kp) :: lnRhoReh,chi,eps1,eps2,eps3

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(3) :: vecbuffer

  logical, parameter :: display = .true.
  logical, parameter :: inversion = .true.

  p = 2._kp
  mu = 1_kp
  wreh = -0.3_kp

  kksfData%real1 = p
  kksfData%real2 = mu
  kksfData%real3 = wreh

  Pstar = powerAmpScalar


  call delete_file('kksfreh_eps2eps1.dat')

!chi stands for phi/mu
 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = kksf_lnrhoend(p,mu,Pstar)

  print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     chi = kksf_x_reheat(p,mu,wreh,lnRhoReh,Pstar,bfold)     
     eps1 = kksf_epsilon_one(chi,p,mu)
     eps2 = kksf_epsilon_two(chi,p,mu)

     if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfold),eps1

     call livewrite('kksfreh_eps2eps1.dat',eps2,eps1,abs(bfold),lnRhoReh)
  
  end do

 
end program kksftest
