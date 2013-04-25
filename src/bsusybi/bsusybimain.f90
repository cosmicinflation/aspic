!test the reheating derivation from slow-roll
program bsusybimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use bsusybisr, only : bsusybi_epsilon_one, bsusybi_epsilon_two 
  use bsusybisr, only : bsusybi_epsilon_three,bsusybi_xendmax
  use bsusybireheat, only : bsusybi_lnrhoreh_max, bsusybi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use bsusybisr, only : bsusybi_norm_potential
  use bsusybireheat, only : bsusybi_x_rreh, bsusybi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,ngamma,nxend

  real(kp) :: gammaBSUSYB,xend,w,bfoldstar,gammamin,gammamax,xendmin,xendmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End

  Pstar = powerAmpScalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!          Calculates the prior space and               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ngamma=1000
  gammamin=10._kp**(-5._kp)
  gammamax=1._kp/sqrt(3._kp)*0.99_kp

  call delete_file('bsusybi_xendmax.dat')
  do i=1,ngamma
     gammaBSUSYB=gammamin+(gammamax-gammamin)*(real(i,kp)/real(ngamma,kp))

     call livewrite('bsusybi_xendmax.dat',gammaBSUSYB,bsusybi_xendmax(40._kp,gammaBSUSYB), &
          bsusybi_xendmax(60._kp,gammaBSUSYB),bsusybi_xendmax(80._kp,gammaBSUSYB))
  end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts = 20
  ngamma=10
  nxend=10

  gammamin=10._kp**(-3._kp)
  gammamax=1._kp/sqrt(3._kp)*0.5_kp

  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('bsusybi_predic.dat')
  call delete_file('bsusybi_nsr.dat')


  do j=1,ngamma
     gammaBSUSYB=exp(log(gammamin) +(log(gammamax)-log(gammamin)) &
          * real(j-1,kp)/real(ngamma-1,kp))

     !Prior on xend

     xendmax=bsusybi_xendmax(70._kp,gammaBSUSYB)
     xendmin=2._kp*xendmax


     do k=1,nxend
        xend=xendmin*(xendmax/xendmin)**(real(k,kp)/real(nxend,kp))

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = bsusybi_lnrhoreh_max(gammaBSUSYB,xend,Pstar)

        print *,'gamma=',gammaBSUSYB,'xend=',xend,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = bsusybi_x_star(gammaBSUSYB,xend,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = bsusybi_epsilon_one(xstar,gammaBSUSYB)
           eps2 = bsusybi_epsilon_two(xstar,gammaBSUSYB)
           eps3 = bsusybi_epsilon_three(xstar,gammaBSUSYB)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('bsusybi_predic.dat',gammaBSUSYB,xend,xend/xendmax,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('bsusybi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  gammaBSUSYB = 1e-2
  xend = bsusybi_xendmax(70._kp,gammaBSUSYB)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = bsusybi_x_rrad(gammaBSUSYB,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = bsusybi_epsilon_one(xstar,gammaBSUSYB)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  bsusybi_epsilon_one(xend,gammaBSUSYB)
     VendOverVstar = bsusybi_norm_potential(xend,gammaBSUSYB) &
          /bsusybi_norm_potential(xstar,gammaBSUSYB)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = bsusybi_x_rreh(gammaBSUSYB,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = bsusybi_x_star(gammaBSUSYB,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program bsusybimain
