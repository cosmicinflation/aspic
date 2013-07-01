!test the reheating derivation from slow-roll
program cncimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use cncisr, only : cnci_epsilon_one, cnci_epsilon_two,cnci_epsilon_three
  use cncisr, only : cnci_xendmin, cnci_x_epsoneunity
  use cncireheat, only : cnci_lnrhoreh_max, cnci_x_star
  use infinout, only : delete_file, livewrite, has_not_shifted
  use srreheat, only : log_energy_reheat_ingev

  use cncisr, only : cnci_norm_potential
  use cncireheat, only : cnci_x_rreh, cnci_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat


  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nxend

  real(kp) :: alpha,xend,w,bfoldstar,alphamin,alphamax,xendmin,xendmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r
  real(kp), dimension(10) ::alphavalues

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End

  Pstar = powerAmpScalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            Calculates the prior space                 !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nalpha=1000
  alphamin=10._kp**(-5._kp)
  alphamax=1._kp

  call delete_file('cnci_xendmax.dat')
  do i=1,nalpha
     alpha=alphamin+(alphamax-alphamin)*(real(i,kp)/real(nalpha,kp))

     call livewrite('cnci_xendmin.dat',alpha,cnci_xendmin(40._kp,alpha), &
          cnci_xendmin(60._kp,alpha),cnci_xendmin(80._kp,alpha))
  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts = 20
  nxend=20


  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('cnci_predic.dat')
  call delete_file('cnci_nsr.dat')

  alphavalues(3)=10._kp**(-3.)
  alphavalues(2)=10._kp**(-1.)
  alphavalues(1)=0.2_kp

  do j=1,3
     alpha=alphavalues(j)



     !Prior on xend
     if (alpha .eq. 10._kp**(-3.)) then
        xendmin=cnci_xendmin(55._kp,alpha)
        xendmax=30._kp*xendmin
        nxend=400
     endif
     if (alpha .eq. 10._kp**(-1.)) then
        xendmin=cnci_xendmin(58._kp,alpha)
        xendmax=10._kp*xendmin
        nxend=100
     endif
     if (alpha .eq. 0.2_kp) then
        xendmin=cnci_xendmin(60._kp,alpha)
        xendmax=2._kp*xendmin
        nxend=100
     endif



     do k=1,nxend
        xend=xendmin+(xendmax-xendmin)*(real(k,kp)/real(nxend,kp))

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = cnci_lnrhoreh_max(alpha,xend,Pstar)

        print *,'alpha=',alpha,'xend=',xend,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = cnci_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = cnci_epsilon_one(xstar,alpha)
           eps2 = cnci_epsilon_two(xstar,alpha)
           eps3 = cnci_epsilon_three(xstar,alpha)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           if (has_not_shifted(0.0075_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
              cycle
           endif
          
           call livewrite('cnci_predic.dat',alpha,xend,xendmin,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('cnci_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 1e-2
  xend = 100
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = cnci_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = cnci_epsilon_one(xstar,alpha)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  cnci_epsilon_one(xend,alpha)
     VendOverVstar = cnci_norm_potential(xend,alpha)/cnci_norm_potential(xstar,alpha)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = cnci_x_rreh(alpha,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = cnci_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program cncimain
