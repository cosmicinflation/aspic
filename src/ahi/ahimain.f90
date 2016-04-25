!test the reheating derivation from slow-roll
program ahimain
  use infprec, only : kp, pi
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ahisr, only : ahi_epsilon_one, ahi_epsilon_two, ahi_epsilon_three
  use ahireheat, only : ahi_lnrhoreh_max, ahi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use ahisr, only : ahi_norm_potential, ahi_x_endinf
  use ahireheat, only : ahi_x_rreh, ahi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat


  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  real(kp) :: phi0,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(1:4) ::phi0values

  real(kp)  :: alpha,alphamin,alphamax,eps1A,eps2A,eps3A,nsA,rA
  real(kp)  :: eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  integer :: nalpha

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  phi0values(1)=10._kp**(0.)
  phi0values(2)=3.*10._kp**(0.)
  phi0values(3)=10._kp**(1.)
  phi0values(4)=10._kp**(2.)

  Pstar = powerAmpScalar


  call delete_file('ahi_predic.dat')
  call delete_file('ahi_nsr.dat')

  do j=1,size(phi0values)

     phi0=phi0values(j)

     w=0._kp

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = ahi_lnrhoreh_max(phi0,Pstar)

     print *,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = ahi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = ahi_epsilon_one(xstar,phi0)
        eps2 = ahi_epsilon_two(xstar,phi0)
        eps3 = ahi_epsilon_three(xstar,phi0)


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(pi**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('ahi_predic.dat',phi0,eps1,eps2,eps3,r,ns,Treh)

     end do
  end do


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  phi0 = 1.
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = ahi_x_rrad(phi0,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = ahi_epsilon_one(xstar,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = ahi_x_endinf(phi0)
     eps1end =  ahi_epsilon_one(xend,phi0)
     VendOverVstar = ahi_norm_potential(xend,phi0)/ahi_norm_potential(xstar,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = ahi_x_rreh(phi0,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = ahi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program ahimain
