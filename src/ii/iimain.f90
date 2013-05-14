!test the reheating derivation from slow-roll
program iimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use iisr, only : ii_epsilon_one, ii_epsilon_two, ii_epsilon_three, ii_xendmin
  use iireheat, only : ii_lnrhoreh_max, ii_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use iisr, only : ii_norm_potential
  use iireheat, only : ii_x_rreh, ii_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts

  real(kp) :: beta,xendinf,w,bfoldstar,xendmin,xendmax,betamin,betamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar

  !Calculates the prior space
  npts=1000
  betamin=epsilon(1._kp)
  betamax=60._kp

  call delete_file('ii_prior.dat')

  do i=1,npts
     beta=betamin+(betamax-betamin)*(real(i,kp)/real(npts,kp))
     call livewrite('ii_prior.dat',beta,ii_xendmin(60._kp,beta))
  end do

  !Calculates the CMB constraints space

  npts = 20

  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('ii_predic.dat')
  call delete_file('ii_nsr.dat')

  betamin=0.1_kp
  betamax=10._kp

  betamin=1._kp
  betamax=70._kp

  do j=0,3
     beta=betamin*(betamax/betamin)**(real(j,kp)/real(3,kp))

     xendmin=ii_xendmin(60._kp,beta)
     xendmax=100._kp*xendmin

     do k=1,20
        xendinf=xendmin*(xendmax/xendmin)**(real(k,kp)/real(20,kp))

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = ii_lnrhoreh_max(beta,xendinf,Pstar)

        print *,'beta=',beta,'xendinf=',xendinf,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = ii_x_star(beta,xendinf,w,lnRhoReh,Pstar,bfoldstar)

           !print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

           eps1 = ii_epsilon_one(xstar,beta)
           eps2 = ii_epsilon_two(xstar,beta)
           eps3 = ii_epsilon_three(xstar,beta)

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('ii_predic.dat',beta,xendinf,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('ii_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  beta = 0.5
  xend = 20
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = ii_x_rrad(beta,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = ii_epsilon_one(xstar,beta)

!consistency test
!get lnR from lnRrad and check that it gives the same xstar
     eps1end =  ii_epsilon_one(xend,beta)
     VendOverVstar = ii_norm_potential(xend,beta)/ii_norm_potential(xstar,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = ii_x_rreh(beta,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

!second consistency check
!get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = ii_x_star(beta,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program iimain
