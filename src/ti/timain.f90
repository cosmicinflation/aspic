!test the reheating derivation from slow-roll
program timain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use tisr, only : ti_epsilon_one, ti_epsilon_two,ti_epsilon_three
  use tireheat, only : ti_lnrhoreh_max, ti_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use tisr, only : ti_norm_potential, ti_x_endinf
  use tireheat, only : ti_x_rreh, ti_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nmu

  real(kp) :: alpha,mu,w,bfoldstar,alphamin,alphamax,mumin,mumax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(:), allocatable ::muvalues

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(muvalues(1:3))

  muvalues(1)=10._kp**(-6.)
  muvalues(2)=10._kp**(-4.)
  muvalues(3)=10._kp**(-2.)

  npts = 6
  nalpha=40

  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('ti_predic.dat')
  call delete_file('ti_nsr.dat')

  do j=1,size(muvalues)
     mu=muvalues(j)

     alphamin=0.5_kp-mu**2/10._kp
     alphamax=0.5_kp+mu**2/20._kp

     do k=0,nalpha
        alpha=alphamin+(alphamax-alphamin)*(real(k,kp)/real(nalpha,kp))


        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = ti_lnrhoreh_max(alpha,mu,Pstar)

        print *,'alpha=',alpha,'mu/Mp=',mu,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = ti_x_star(alpha,mu,w,lnRhoReh,Pstar,bfoldstar)

           eps1 = ti_epsilon_one(xstar,alpha,mu)
           eps2 = ti_epsilon_two(xstar,alpha,mu)
           eps3 = ti_epsilon_three(xstar,alpha,mu)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('ti_predic.dat',alpha,(1-2._kp*alpha)/mu**2,mu,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('ti_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do

  nmu=80
  mumin=10._kp**(-7.)
  mumax=10._kp**(0.)
  alpha=0.5_kp

  do j=1,nmu
     mu=mumin*(mumax/mumin)**(real(j,kp)/real(nmu,kp))

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = ti_lnrhoreh_max(alpha,mu,Pstar)

     print *,'alpha=',alpha,'mu/Mp=',mu,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = ti_x_star(alpha,mu,w,lnRhoReh,Pstar,bfoldstar)

        eps1 = ti_epsilon_one(xstar,alpha,mu)
        eps2 = ti_epsilon_two(xstar,alpha,mu)
        eps3 = ti_epsilon_three(xstar,alpha,mu)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('ti_predic.dat',alpha,0._kp,mu,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('ti_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 0.5
  mu = 0.9
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = ti_x_rrad(alpha,mu,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = ti_epsilon_one(xstar,alpha,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = ti_x_endinf(alpha,mu)
     eps1end =  ti_epsilon_one(xend,alpha,mu)
     VendOverVstar = ti_norm_potential(xend,alpha,mu)/ti_norm_potential(xstar,alpha,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = ti_x_rreh(alpha,mu,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = ti_x_star(alpha,mu,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo




end program timain
