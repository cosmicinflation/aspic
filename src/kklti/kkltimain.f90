!test the reheating derivation from slow-roll
program bimain
  use infprec, only : kp, transfert
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use kkltireheat, only : kklti_x_star, kklti_lnrhoreh_max
  use kkltisr, only : kklti_epsilon_one, kklti_epsilon_two, kklti_epsilon_three
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use kkltisr, only : kklti_norm_potential, kklti_x_endinf
  use kkltireheat, only : kklti_x_rreh, kklti_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none

  real(kp) :: Pstar

  integer :: i,j
  integer :: npts = 20,nmu=30.

  real(kp) :: mu,p,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r,Treh
  real(kp) :: logErehGeV

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp) ::mumin,mumax

  logical, parameter :: display = .true.

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  w = 0._kp

  Pstar = powerAmpScalar

  call delete_file('kklti_predic.dat')
  call delete_file('kklti_nsr.dat')

!!!!!!!!!!!!!! 
!!!! p=2  !!!!
!!!!!!!!!!!!!!

  mumin=10._kp**(-3.)
  mumax=10._kp**(3.)
  p = 2._kp


  do j=0,nmu

     ! Logarithmic step
     mu=mumin*(mumax/mumin)**(real(j,kp)/real(nmu,kp))

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = kklti_lnrhoreh_max(p,mu,Pstar)

     print *,'lnRhoRehMin= ',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = kklti_x_star(p,mu,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = kklti_epsilon_one(xstar,p,mu)
        eps2 = kklti_epsilon_two(xstar,p,mu)
        eps3 = kklti_epsilon_three(xstar,p,mu)

        if (display) print *,'lnRhoReh=',lnRhoReh,'  N*=',abs(bfoldstar),'eps1star=',eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('kklti_predic.dat',p,mu,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('kklti_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do

!!!!!!!!!!!!!! 
!!!! p=3  !!!!
!!!!!!!!!!!!!!

  mumin=10._kp**(-3.)
  mumax=10._kp**(3.)
  p = 3._kp


  do j=0,nmu

     ! Logarithmic step
     mu=mumin*(mumax/mumin)**(real(j,kp)/real(nmu,kp))

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = kklti_lnrhoreh_max(p,mu,Pstar)

     print *,'lnRhoRehMin= ',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = kklti_x_star(p,mu,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = kklti_epsilon_one(xstar,p,mu)
        eps2 = kklti_epsilon_two(xstar,p,mu)
        eps3 = kklti_epsilon_three(xstar,p,mu)

        if (display) print *,'lnRhoReh=',lnRhoReh,'  N*=',abs(bfoldstar),'eps1star=',eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('kklti_predic.dat',p,mu,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('kklti_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do

!!!!!!!!!!!!!! 
!!!! p=4  !!!!
!!!!!!!!!!!!!!

  mumin=10._kp**(-3.)
  mumax=10._kp**(3.)
  p = 4._kp


  do j=0,nmu

     ! Logarithmic step
     mu=mumin*(mumax/mumin)**(real(j,kp)/real(nmu,kp))

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = kklti_lnrhoreh_max(p,mu,Pstar)

     print *,'lnRhoRehMin= ',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = kklti_x_star(p,mu,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = kklti_epsilon_one(xstar,p,mu)
        eps2 = kklti_epsilon_two(xstar,p,mu)
        eps3 = kklti_epsilon_three(xstar,p,mu)

        if (display) print *,'lnRhoReh=',lnRhoReh,'  N*=',abs(bfoldstar),'eps1star=',eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('kklti_predic.dat',p,mu,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('kklti_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

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

     xstar = kklti_x_rrad(p,mu,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = kklti_epsilon_one(xstar,p,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = kklti_x_endinf(p,mu)
     eps1end =  kklti_epsilon_one(xend,p,mu)
     VendOverVstar = kklti_norm_potential(xend,p,mu)/kklti_norm_potential(xstar,p,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = kklti_x_rreh(p,mu,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = kklti_x_star(p,mu,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program bimain
