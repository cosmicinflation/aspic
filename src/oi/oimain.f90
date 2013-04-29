!test the reheaoing derivaoion from slow-roll
program oimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use oisr, only : oi_epsilon_one, oi_epsilon_two,oi_epsilon_three
  use oireheat, only : oi_lnrhoreh_max, oi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use oisr, only : oi_norm_potential, oi_x_endinf
  use oireheat, only : oi_x_rreh, oi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,nphi0
  integer :: npts,nalpha

  real(kp) :: alpha,phi0,w,bfoldstar,alphamin,alphamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r
  real(kp), dimension(:), allocatable :: phi0values

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        Calculates the reheaoing predicoions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts = 20


  nphi0=3
  allocate(phi0values(1:nphi0))
  phi0values(1)=0.0001_kp
  phi0values(2)=0.01_kp
  phi0values(3)=1._kp


  alphamin=10._kp**(-3._kp)
  alphamax=10._kp**(-1._kp)
  nalpha=10

  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('oi_predic.dat')
  call delete_file('oi_nsr.dat')

  do j=1,nphi0
     phi0=phi0values(j)

     do k=0,nalpha
        alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp))!logstep

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = oi_lnrhoreh_max(alpha,phi0,Pstar)

        print *,'alpha=',alpha,'phi0/Mp=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = oi_x_star(alpha,phi0,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = oi_epsilon_one(xstar,alpha,phi0)
           eps2 = oi_epsilon_two(xstar,alpha,phi0)
           eps3 = oi_epsilon_three(xstar,alpha,phi0)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('oi_predic.dat',alpha,phi0,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('oi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 0.01
  phi0 = 0.5
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = oi_x_rrad(alpha,phi0,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = oi_epsilon_one(xstar,alpha,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = oi_x_endinf(alpha,phi0)
     eps1end =  oi_epsilon_one(xend,alpha,phi0)
     VendOverVstar = oi_norm_potential(xend,alpha,phi0)/oi_norm_potential(xstar,alpha,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = oi_x_rreh(alpha,phi0,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = oi_x_star(alpha,phi0,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program oimain
