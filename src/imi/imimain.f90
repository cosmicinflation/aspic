!test the reheating derivation from slow-roll
program imimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use imisr, only : imi_epsilon_one, imi_epsilon_two, imi_epsilon_three, imi_xendmin
  use imireheat, only : imi_lnrhoreh_max, imi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use imisr, only : imi_norm_potential
  use imireheat, only : imi_x_rreh, imi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,nxend
  integer :: npts = 8

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: pmin,pmax,p,xendmin,xendmax,xend

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End

  Pstar = powerAmpScalar

  nxend=12

  w=0._kp

  call delete_file('imi_predic.dat')
  call delete_file('imi_nsr.dat')

  pmin=1.
  pmax=6.

  do j=0,int(pmax-pmin)

     p = pmin+real(j,kp)

     xendmin=imi_xendmin(65._kp,p)
     xendmax=12._kp*xendmin

     do k=0,nxend

        xend=xendmin*(xendmax/xendmin)**(real(k,kp)/real(nxend,kp))

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = imi_lnrhoreh_max(p,xend,Pstar)

        print *,'p=',p,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = imi_x_star(p,xend,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar

           eps1 = imi_epsilon_one(xstar,p)
           eps2 = imi_epsilon_two(xstar,p)
           eps3 = imi_epsilon_three(xstar,p)

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('imi_predic.dat',p,xend,eps1,eps2,eps3,r,ns,Treh)


        end do

     end do

  end do


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  p = 4
  xend = 20*imi_xendmin(65._kp,p)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = imi_x_rrad(p,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = imi_epsilon_one(xstar,p)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  imi_epsilon_one(xend,p)
     VendOverVstar = imi_norm_potential(xend,p)/imi_norm_potential(xstar,p)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = imi_x_rreh(p,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = imi_x_star(p,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program imimain
