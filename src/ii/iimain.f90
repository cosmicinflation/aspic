!test the reheating derivation from slow-roll
program iimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use iisr, only : ii_epsilon_one, ii_epsilon_two, ii_epsilon_three, ii_prior_xendmin
  use iireheat, only : ii_lnrhoend, ii_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts

  real(kp) :: beta,xendinf,w,bfoldstar,xendmin,xendmax,betamin,betamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer


  Pstar = powerAmpScalar

!Calculates the prior space
  npts=1000
  betamin=epsilon(1._kp)
  betamax=20._kp

  call delete_file('ii_prior.dat')

  do i=1,npts
       beta=betamin*(betamax/betamin)**(real(i,kp)/real(npts,kp))
       call livewrite('ii_prior.dat',beta,ii_prior_xendmin(beta,-60._kp))
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

    xendmin=ii_prior_xendmin(beta,-60._kp)
    xendmax=100._kp*xendmin

    do k=1,20
       xendinf=xendmin*(xendmax/xendmin)**(real(k,kp)/real(20,kp))

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = ii_lnrhoend(beta,xendinf,Pstar)

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



end program iimain
