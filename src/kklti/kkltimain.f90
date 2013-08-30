!test the reheating derivation from slow-roll
program kkltimain
  use infprec, only : pi, kp, transfert
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use kkltireheat, only : kklti_x_star, kklti_lnrhoreh_max
  use kkltisr, only : kklti_epsilon_one, kklti_epsilon_two, kklti_epsilon_three
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use kkltisr, only : kklti_norm_potential, kklti_x_epsoneunity
  use kkltisr, only : kklti_ln_xstg, kklti_ln_xuv, kklti_x_trajectory
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

  real(kp) :: alpha, calN, v, gstring, vbar, y
  real(kp) :: xstg, xuv, xeps, xini, tachEnd

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

     xend = kklti_x_epsoneunity(p,mu)

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = kklti_lnrhoreh_max(p,mu,xend,Pstar)

     print *,'lnRhoRehMin= ',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = kklti_x_star(p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
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

     xend = kklti_x_epsoneunity(p,mu)

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = kklti_lnrhoreh_max(p,mu,xend,Pstar)

     print *,'lnRhoRehMin= ',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = kklti_x_star(p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
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

     xend = kklti_x_epsoneunity(p,mu)

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = kklti_lnrhoreh_max(p,mu,xend,Pstar)

     print *,'lnRhoRehMin= ',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = kklti_x_star(p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
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


!!!!!!!!!!!!!!!!!!!! 
!!!! p=4 STRING !!!!
!!!!!!!!!!!!!!!!!!!!

  call delete_file('kkltistg_nsr.dat')
  call delete_file('kkltistg_predic.dat')

  mumin=10._kp**(-6._kp)
  mumax=0.1_kp
  p = 4._kp
  
  calN = 5
  v =16._kp/27_kp
  gstring = 0.005_kp
  alpha = 0.25_kp

  y = 4._kp*pi*gstring*calN/v
  vbar = v/(4._kp*pi*gstring)**2
  


  mumin=10._kp**(-6._kp)
  mumax=0.1_kp
  p = 4._kp


  do j=0,nmu

     ! Logarithmic step
     mu=mumin*(mumax/mumin)**(real(j,kp)/real(nmu,kp))

     xeps = kklti_x_epsoneunity(p,mu)
     xstg = exp(kklti_ln_xstg(p,mu,y,vbar,calN,alpha))
     xuv = exp(kklti_ln_xuv(p,mu,y,vbar,calN,alpha))

     print *,'string BI xeps xstg xuv',xeps,xstg,xuv

     xend = max(xeps,xstg)
     if (xeps.gt.xstg) then
       tachEnd=1._kp
     else
       tachEnd=0._kp
     endif
     xini = kklti_x_trajectory(-110._kp,xend,p,mu)

     
     if (xini.gt.xuv) then
        print *,'UV bound violated: skipping xini xuv',xini,xuv
        cycle
     endif

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = kklti_lnrhoreh_max(p,mu,xend,Pstar)

     print *,'lnRhoRehMin= ',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = kklti_x_star(p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = kklti_epsilon_one(xstar,p,mu)
        eps2 = kklti_epsilon_two(xstar,p,mu)
        eps3 = kklti_epsilon_three(xstar,p,mu)

        if (display) print *,'lnRhoReh=',lnRhoReh,'  N*=',abs(bfoldstar),'eps1star=',eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('kkltistg_predic.dat',p,mu,eps1,eps2,eps3,r,ns,Treh,tachEnd)

        call livewrite('kkltistg_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

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
     
     xend = kklti_x_epsoneunity(p,mu)
     xstar = kklti_x_rrad(p,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = kklti_epsilon_one(xstar,p,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar

     eps1end =  kklti_epsilon_one(xend,p,mu)
     VendOverVstar = kklti_norm_potential(xend,p,mu)/kklti_norm_potential(xstar,p,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = kklti_x_rreh(p,mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = kklti_x_star(p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program kkltimain
