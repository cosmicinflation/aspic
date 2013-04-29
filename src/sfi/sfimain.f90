!test the reheating derivation from slow-roll
program sfimain
  use infprec, only : kp, transfert
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use sfireheat, only : sfi_x_star, sfi_lnrhoreh_max
  use sfireheat, only : sfi_xpmu_fromepsilon, sfi_lnrhoreh_fromepsilon
  use sfisr, only : sfi_epsilon_one, sfi_epsilon_two,sfi_epsilon_three
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use sfisr, only : sfi_norm_potential, sfi_x_endinf
  use sfireheat, only : sfi_x_rreh, sfi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none

  type(transfert) :: sfiData
  real(kp) :: Pstar,calF

  integer :: i,j
  integer :: npts = 20,nmu=30.

  real(kp) :: mu,p,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r,Treh
  real(kp) :: logErehGeV

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(3) :: vecbuffer

  real(kp) ::mumin,mumax

  logical, parameter :: display = .true.
  logical, parameter :: inversion = .true.

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  w = 0._kp

  Pstar = powerAmpScalar

  call delete_file('sfi_predic.dat')
  call delete_file('sfi_nsr.dat')

!!!!!!!!!!!!!! 
!!!! p=1  !!!!
!!!!!!!!!!!!!!

  mumin=12.
  mumax=10._kp**(3.)
  p = 1._kp


  do j=0,nmu

     ! Ultralogarithmic step
     mu=mumin/exp(1._kp)*exp(exp(real(j,kp)/real(nmu,kp)*log(1._kp+log(mumax/mumin))))


     !xstar stands for phistar/mu

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = sfi_lnrhoreh_max(p,mu,Pstar)

     print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = sfi_x_star(p,mu,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = sfi_epsilon_one(xstar,p,mu)
        eps2 = sfi_epsilon_two(xstar,p,mu)
        eps3 = sfi_epsilon_three(xstar,p,mu)

        if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfoldstar),eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('sfi_predic.dat',p,mu,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('sfi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do

!!!!!!!!!!!!!! 
!!!! p=1  !!!!
!!!!!!!!!!!!!!

  mumin=0.7
  mumax=10._kp**(3.)
  p = 2._kp

  do j=0,nmu


     ! Ultralogarithmic step
     mu=mumin/exp(1._kp)*exp(exp(real(j,kp)/real(nmu,kp)*log(1._kp+log(mumax/mumin))))


     !xstar stands for phistar/mu

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = sfi_lnrhoreh_max(p,mu,Pstar)

     print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = sfi_x_star(p,mu,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = sfi_epsilon_one(xstar,p,mu)
        eps2 = sfi_epsilon_two(xstar,p,mu)
        eps3 = sfi_epsilon_three(xstar,p,mu)

        if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfoldstar),eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('sfi_predic.dat',p,mu,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('sfi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do

!!!!!!!!!!!!!! 
!!!! p=1  !!!!
!!!!!!!!!!!!!!

  mumin=0.7
  mumax=10._kp**(3.)
  p = 4._kp

  do j=0,nmu


     ! Ultralogarithmic step
     mu=mumin/exp(1._kp)*exp(exp(real(j,kp)/real(nmu,kp)*log(1._kp+log(mumax/mumin))))


     !xstar stands for phistar/mu

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = sfi_lnrhoreh_max(p,mu,Pstar)

     print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = sfi_x_star(p,mu,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = sfi_epsilon_one(xstar,p,mu)
        eps2 = sfi_epsilon_two(xstar,p,mu)
        eps3 = sfi_epsilon_three(xstar,p,mu)

        if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfoldstar),eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('sfi_predic.dat',p,mu,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('sfi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  p=2.5
  mu = 0.1
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = sfi_x_rrad(p,mu,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = sfi_epsilon_one(xstar,p,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = sfi_x_endinf(p,mu)
     eps1end =  sfi_epsilon_one(xend,p,mu)
     VendOverVstar = sfi_norm_potential(xend,p,mu)/sfi_norm_potential(xstar,p,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = sfi_x_rreh(p,mu,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = sfi_x_star(p,mu,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program sfimain
