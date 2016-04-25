!test the reheating derivation from slow-roll
program sdimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use sdisr, only : sdi_epsilon_one, sdi_epsilon_two, sdi_epsilon_three, sdi_check_params
  use sdireheat, only : sdi_lnrhoreh_max, sdi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use sdisr, only : sdi_norm_potential
  use sdireheat, only : sdi_x_rreh, sdi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts, Nphi0, Nxend

  real(kp) :: phi0,xendinf,w,bfoldstar,xendmin,xendmax,phi0min,phi0max
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar




  npts = 20

  w=0._kp
  !  w = 1._kp/3._kp


  call delete_file('sdi_predic.dat')
  call delete_file('sdi_nsr.dat')

  phi0min=1._kp/sqrt(2._kp)
  phi0min=5._kp
  phi0max=10._kp**3

  Nphi0 = 30

  xendmin=10._kp**(-3.)
  xendmax=50._kp

  Nxend = 200


  do j=0,Nphi0
     phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(Nphi0,kp))

     do k=1,Nxend
        xendinf=xendmin*(xendmax/xendmin)**(real(k,kp)/real(Nxend,kp))

        if (sdi_check_params(phi0, xendinf)) then

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = sdi_lnrhoreh_max(phi0,xendinf,Pstar)

        print *,'phi0=',phi0,'xendinf=',xendinf,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = sdi_x_star(phi0,xendinf,w,lnRhoReh,Pstar,bfoldstar)

           !print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

           eps1 = sdi_epsilon_one(xstar,phi0)
           eps2 = sdi_epsilon_two(xstar,phi0)
           eps3 = sdi_epsilon_three(xstar,phi0)

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           print*, 'ns=',ns,'r=',r,'bfoldstar=',bfoldstar

           call livewrite('sdi_predic.dat',phi0,xendinf,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('sdi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

      end if

     end do

  end do

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  phi0 = 100._kp
  xend = 20._kp
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = sdi_x_rrad(phi0,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = sdi_epsilon_one(xstar,phi0)

!consistency test
!get lnR from lnRrad and check that it gives the same xstar
     eps1end =  sdi_epsilon_one(xend,phi0)
     VendOverVstar = sdi_norm_potential(xend)/sdi_norm_potential(xstar)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = sdi_x_rreh(phi0,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

!second consistency check
!get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = sdi_x_star(phi0,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program sdimain
