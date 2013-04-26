!test the reheating derivation from slow-roll
program kmiimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use kmiisr, only : kmii_epsilon_one, kmii_epsilon_two, kmii_epsilon_three
  use kmiireheat, only : kmii_lnrhoreh_max, kmii_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use kmiisr, only : kmii_norm_potential, kmii_x_endinf
  use kmiireheat, only : kmii_x_rreh, kmii_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer
  real(kp) ::alphamin,alphamax

  real(kp) :: eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  integer :: nalpha

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  alphamin=sqrt(2.)/(sqrt(2.)+1.)*exp((2.+sqrt(2.))/(1.+sqrt(2.)))
  alphamax=exp(1.)
  alphamax=10000._kp

  Pstar = powerAmpScalar
  !  w = 1._kp/3._kp
  w=0._kp

  nalpha=10

  call delete_file('kmii_predic.dat')
  call delete_file('kmii_nsr.dat')

  do j=0,nalpha+1
     alpha = alphamin+(real(j-1,kp)/real(nalpha,kp))*(alphamax-alphamin) !arithmetic step
     alpha = alphamin*(alphamax/alphamin)**(real(j-1,kp)/real(nalpha,kp)) !logarithmic step


     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = kmii_lnrhoreh_max(alpha,Pstar)

     print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = kmii_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = kmii_epsilon_one(xstar,alpha)
        eps2 = kmii_epsilon_two(xstar,alpha)
        eps3 = kmii_epsilon_three(xstar,alpha)

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('kmii_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('kmii_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('kmii_predic_summarized.dat') 
  nalpha=1000
  alphamin=sqrt(2.)/(sqrt(2.)+1.)*exp((2.+sqrt(2.))/(1.+sqrt(2.)))
  alphamax=exp(1.)
  w=0._kp
  do j=1,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))
     lnRhoReh = lnRhoNuc
     xstarA = kmii_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     eps1A = kmii_epsilon_one(xstarA,alpha)
     eps2A = kmii_epsilon_two(xstarA,alpha)
     eps3A = kmii_epsilon_three(xstarA,alpha)
     nsA = 1._kp - 2._kp*eps1A - eps2A
     rA = 16._kp*eps1A
     lnRhoReh = kmii_lnrhoreh_max(alpha,Pstar)
     xstarB = kmii_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     eps1B = kmii_epsilon_one(xstarB,alpha)
     eps2B = kmii_epsilon_two(xstarB,alpha)
     eps3B = kmii_epsilon_three(xstarB,alpha)
     nsB = 1._kp - 2._kp*eps1B - eps2B
     rB =16._kp*eps1B
     call livewrite('kmii_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
  enddo

 write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 10
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = kmii_x_rrad(alpha,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = kmii_epsilon_one(xstar,alpha)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = kmii_x_endinf(alpha)
     eps1end =  kmii_epsilon_one(xend,alpha)
     VendOverVstar = kmii_norm_potential(xend,alpha)/kmii_norm_potential(xstar,alpha)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = kmii_x_rreh(alpha,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = kmii_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program kmiimain
