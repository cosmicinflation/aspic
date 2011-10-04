!test the reheating derivation from slow-roll
program sfrehmain  
  use infprec, only : kp, transfert
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use sfreheat, only : sf_x_reheat, sf_lnrhoend
  use sfreheat, only : sf_x_obs, sf_lnrhoreh
  use sfsrevol, only : sf_epsilon_one, sf_epsilon_two
  use infinout, only : delete_file, livewrite
  implicit none

  type(transfert) :: sfData
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
  mu = 8_kp
  wreh = -0.3_kp

  sfData%real1 = p
  sfData%real2 = mu
  sfData%real3 = wreh

  Pstar = powerAmpScalar

  call delete_file('sfreh_eps2eps1.dat')

!chi stands for phi/mu
 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = sf_lnrhoend(p,mu,Pstar)

  print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     chi = sf_x_reheat(p,mu,wreh,lnRhoReh,Pstar,bfold)
     eps1 = sf_epsilon_one(chi,p,mu)
     eps2 = sf_epsilon_two(chi,p,mu)

     if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfold),eps1

     call livewrite('sfreh_eps2eps1.dat',eps2,eps1,abs(bfold),lnRhoReh)
  
  end do


  if (inversion) then

     eps1 = 1e-3
     eps2 = 0.04
     eps3 = 0.001
     vecbuffer = sf_x_obs(eps1,eps2,eps3,bfold)

     chi = vecbuffer(1)
     p = vecbuffer(2)
     mu = vecbuffer(3)

     print *,'p= mu= ',p,mu
     print *,'chistar bfold', chi,bfold


     lnRhoReh = sf_lnrhoreh(wreh,eps1,eps2,eps3,Pstar)
     print *,'lnRhoReh =',lnrhoReh
  endif

end program sfrehmain
