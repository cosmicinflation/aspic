!test the reheating derivation from slow-roll
program nckimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use nckisr, only : ncki_epsilon_one, ncki_epsilon_two, ncki_epsilon_three
  use nckireheat, only : ncki_lnrhoreh_max, ncki_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use nckisr, only : ncki_norm_potential, ncki_x_endinf
  use nckireheat, only : ncki_x_rreh, ncki_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nbeta

  real(kp) :: alpha,beta,w,bfoldstar,alphamin,alphamax,betamin,betamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('ncki_predic.dat')
  call delete_file('ncki_nsr.dat')


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                 Case beta>0                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  npts = 20
  nalpha=15
  nbeta=7

  do j=0,nbeta
     if (j.eq.0) then
        beta=10._kp**(-4._kp)
        alphamax=10._kp**(0._kp)
        alphamin=10._kp**(-4_kp)
     else if (j.eq.1) then
        beta=10._kp**(-2.5_kp)
        alphamax=10._kp**(0._kp)
        alphamin=10._kp**(-4_kp)
     else if (j.eq.2) then
        beta=10._kp**(-2.2_kp)
        alphamax=10._kp**(0._kp)
        alphamin=10._kp**(-4_kp)
     else if (j.eq.3) then 
        beta=10._kp**(-1.9_kp)
        alphamax=10._kp**(0._kp)
        alphamin=10._kp**(-4_kp)
     else if (j.eq.4) then
        beta=10._kp**(-1.7_kp)
        alphamax=10._kp**(0._kp)
        alphamin=10._kp**(-4.7_kp)
     else if (j.eq.5) then
        beta=10._kp**(-1.5_kp)
        alphamax=10._kp**(0._kp)
        alphamin=10._kp**(-4_kp)
     else if (j.eq.6) then
        nalpha=40
        beta=10._kp**(-1.2_kp)
        alphamax=10._kp**(0._kp)
        alphamin=10._kp**(-4.7_kp)
     else if (j.eq.7) then
        beta=0.5_kp
        alphamin=10._kp**(-3.5_kp)
        alphamax=10._kp**(-0._kp)
     endif

     print*,'beta=',beta

     do k=0,nalpha
        alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !log step


        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = ncki_lnrhoreh_max(alpha,beta,Pstar)

        print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           print*,'lnRhoReh=',lnRhoReh

           xstar = ncki_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)



           eps1 = ncki_epsilon_one(xstar,alpha,beta)
           eps2 = ncki_epsilon_two(xstar,alpha,beta)
           eps3 = ncki_epsilon_three(xstar,alpha,beta)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('ncki_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('ncki_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!                 Case beta<0                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  npts = 20
  nalpha=25
  nbeta=3

  do j=0,nbeta
     if (j.eq.0) then
        beta=-10._kp**(-4._kp)
        alphamax=10._kp**(1._kp)
        alphamin=10._kp**(-3_kp)
     else if (j.eq.1) then
        beta=-10._kp**(-2.5_kp)
        alphamax=10._kp**(1._kp)
        alphamin=10._kp**(-3_kp)
     else if (j.eq.2) then
        beta=-10._kp**(-2.1_kp)
        alphamax=10._kp**(1._kp)
        alphamin=10._kp**(-3_kp)
     else if (j.eq.3) then
        beta=-0.5_kp
        alphamax=10._kp**(3._kp)
        alphamin=10._kp**(1._kp)
     endif

     print*,'beta=',beta


     do k=0,nalpha
        alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !log step


        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = ncki_lnrhoreh_max(alpha,beta,Pstar)

        print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           print*,'lnRhoReh=',lnRhoReh

           xstar = ncki_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)



           eps1 = ncki_epsilon_one(xstar,alpha,beta)
           eps2 = ncki_epsilon_two(xstar,alpha,beta)
           eps3 = ncki_epsilon_three(xstar,alpha,beta)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('ncki_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('ncki_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do


 write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 1.
  beta = 0.5
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = ncki_x_rrad(alpha,beta,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = ncki_epsilon_one(xstar,alpha,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = ncki_x_endinf(alpha,beta)
     eps1end =  ncki_epsilon_one(xend,alpha,beta)
     VendOverVstar = ncki_norm_potential(xend,alpha,beta)/ncki_norm_potential(xstar,alpha,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = ncki_x_rreh(alpha,beta,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = ncki_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program nckimain
