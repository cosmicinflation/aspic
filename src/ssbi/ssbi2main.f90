!test the reheating derivation from slow-roll
program ssbi2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ssbi2sr, only : ssbi2_epsilon_one, ssbi2_epsilon_two, ssbi2_epsilon_three
  use ssbi2reheat, only : ssbi2_lnrhoreh_max, ssbi2_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use ssbi2sr, only : ssbi2_norm_potential, ssbi2_x_endinf
  use ssbi2reheat, only : ssbi2_x_rreh, ssbi2_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 20

  integer :: Nalpha,Nbeta
  real(kp) ::alphamin, alphamax, betamin, betamax, alpha, beta

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!    Slow Roll Predictions   !!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Pstar = powerAmpScalar

  call delete_file('ssbi2_predic.dat')
  call delete_file('ssbi2_nsr.dat')

  Nalpha=100
  Nbeta=100

  !  w = 1._kp/3._kp
  w=0._kp

  betamin=10._kp**(-3._kp)
  betamax=10._kp**(3._kp)

  do j=0,3 

     if (j .eq. 0) then
        beta=-10._kp**(-10._kp) 
     end if
     if (j .eq. 1) then
        beta=-10._kp**(-4._kp) 
     end if
     if (j .eq. 2) then
        beta=-10._kp**(-3._kp) 
     end if
     if (j .eq. 3) then 
        beta=-10._kp**(-2._kp) 
     end if


     alphamin=10._kp**(-5._kp)
     alphamax=10._kp**(-1._kp)

     do k=0,Nalpha 
        alpha=-alphamin*(alphamax/alphamin)**(real(k,kp)/Nalpha)  !logarithmic step

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = ssbi2_lnrhoreh_max(alpha,beta,Pstar)


        print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = ssbi2_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar


           eps1 = ssbi2_epsilon_one(xstar,alpha,beta)
           eps2 = ssbi2_epsilon_two(xstar,alpha,beta)
           eps3 = ssbi2_epsilon_three(xstar,alpha,beta)


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)

           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('ssbi2_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('ssbi2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  beta = -0.01
  alpha = -0.001
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = ssbi2_x_rrad(alpha,beta,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = ssbi2_epsilon_one(xstar,alpha,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = ssbi2_x_endinf(alpha,beta)
     eps1end =  ssbi2_epsilon_one(xend,alpha,beta)
     VendOverVstar = ssbi2_norm_potential(xend,alpha,beta)/ssbi2_norm_potential(xstar,alpha,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = ssbi2_x_rreh(alpha,beta,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = ssbi2_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program ssbi2main
