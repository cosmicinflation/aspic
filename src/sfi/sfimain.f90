!test the reheating derivation from slow-roll
program sfimain
  use infprec, only : kp, transfert
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use sfireheat, only : sfi_x_reheat, sfi_lnrhoend
  use sfireheat, only : sfi_xpmu_fromepsilon, sfi_lnrhoreh_fromepsilon
  use sfisr, only : sfi_epsilon_one, sfi_epsilon_two,sfi_epsilon_three
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  implicit none

  type(transfert) :: sfiData
  real(kp) :: Pstar,calF

  integer :: i
  integer :: npts = 20

  real(kp) :: mu,p,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r
  real(kp) :: logErehGeV

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(3) :: vecbuffer

  logical, parameter :: display = .true.
  logical, parameter :: inversion = .true.

  p = 2._kp
  mu = 8_kp
  w = -0.3_kp

  sfiData%real1 = p
  sfiData%real2 = mu
  sfiData%real3 = w

  Pstar = powerAmpScalar

  call delete_file('sfi_predic.dat')
  call delete_file('sfi_nsr.dat')

!xstar stands for phistar/mu
 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = sfi_lnrhoend(p,mu,Pstar)

  print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = sfi_x_reheat(p,mu,w,lnRhoReh,Pstar,bfoldstar)
     eps1 = sfi_epsilon_one(xstar,p,mu)
     eps2 = sfi_epsilon_two(xstar,p,mu)
     eps3 = sfi_epsilon_three(xstar,p,mu)

     if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfoldstar),eps1


     logErehGeV = log_energy_reheat_ingev(lnRhoReh)

     ns = 1._kp - 2._kp*eps1 - eps2
     r =16._kp*eps1

     call livewrite('sfi_predic.dat',p,mu,eps1,eps2,eps3,r,ns,logErehGeV)

     call livewrite('sfi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

  end do


end program sfimain
