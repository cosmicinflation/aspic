!test the reheating derivation from slow-roll
program lpi2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lpi2sr, only : lpi2_epsilon_one, lpi2_epsilon_two, lpi2_epsilon_three
  use lpi2reheat, only : lpi2_lnrhoreh_max, lpi2_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use lpi2sr, only : lpi2_norm_potential, lpi2_x_endinf
  use lpi2reheat, only : lpi2_x_rreh, lpi2_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nphi0

  real(kp) :: phi0,p,q,w,bfoldstar,phi0min,phi0max,pmin,pmax,qmin,qmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  integer :: npq
  real(kp), dimension(10) :: pvalues, qvalues

  real(kp) ::x,xmin,xmax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('lpi2_predic.dat')
  call delete_file('lpi2_nsr.dat')

  call aspicwrite_header('lpi2',labeps12,labnsr,labbfoldreh,(/'phi0','p   ','q   '/))
  
  w=0._kp
  !  w = 1._kp/3._kp

  npts = 15

  npq=3

  pvalues(1)=1._kp
  qvalues(1)=2._kp

  pvalues(2)=2._kp
  qvalues(2)=2._kp

  pvalues(3)=3._kp
  qvalues(3)=4._kp


  phi0min=10._kp**(0._kp)*5._kp
  phi0max=10._kp**(1.9_kp)

  do k=1,3

     p=pvalues(k)
     q=qvalues(k)

     if (k .eq. 1) nphi0=40
     if (k .eq. 2) nphi0=40
     if (k .eq. 3) nphi0=40


     do j=0,nphi0
        phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(nphi0,kp))

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = lpi2_lnrhoreh_max(p,q,phi0,Pstar)


        print *,'phi0=',phi0,'p=',p,'q=',q,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = lpi2_x_star(p,q,phi0,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = lpi2_epsilon_one(xstar,p,q,phi0)
           eps2 = lpi2_epsilon_two(xstar,p,q,phi0)
           eps3 = lpi2_epsilon_three(xstar,p,q,phi0)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('lpi2_predic.dat',p,q,phi0,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('lpi2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/phi0,p,q/))

        end do

     end do

  end do

  call aspicwrite_end()

 write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  p=3.
  q=4.
  phi0=400

  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = lpi2_x_rrad(p,q,phi0,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = lpi2_epsilon_one(xstar,p,q,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = lpi2_x_endinf(p,q,phi0)
     eps1end =  lpi2_epsilon_one(xend,p,q,phi0)
     VendOverVstar = lpi2_norm_potential(xend,p,q,phi0)/lpi2_norm_potential(xstar,p,q,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = lpi2_x_rreh(p,q,phi0,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = lpi2_x_star(p,q,phi0,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program lpi2main
