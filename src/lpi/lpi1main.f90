!test the reheating derivation from slow-roll
program lpi1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lpi1sr, only : lpi1_epsilon_one, lpi1_epsilon_two, lpi1_epsilon_three
  use lpi1reheat, only : lpi1_lnrhoreh_max, lpi1_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use lpi1sr, only : lpi1_norm_potential, lpi1_x_endinf
  use lpi1reheat, only : lpi1_x_rreh, lpi1_x_rrad
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

  call delete_file('lpi1_predic.dat')
  call delete_file('lpi1_nsr.dat')

  call aspicwrite_header('lpi1',labeps12,labnsr,labbfoldreh,(/'phi0','q   ','p   '/))
  
  w=0._kp
  !  w = 1._kp/3._kp

  npts = 20

  npq=3

  pvalues(1)=4._kp
  qvalues(1)=1._kp

  pvalues(2)=4._kp
  qvalues(2)=2._kp

  pvalues(3)=4._kp
  qvalues(3)=3._kp


  phi0min=10._kp**(-3._kp)
  phi0max=10._kp**(3._kp)

  do k=1,3

     p=pvalues(k)
     q=qvalues(k)

     if (k .eq. 1) nphi0=100
     if (k .eq. 2) nphi0=100
     if (k .eq. 3) nphi0=100


     do j=0,nphi0
        phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(nphi0,kp))

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = lpi1_lnrhoreh_max(p,q,phi0,Pstar)


        print *,'phi0=',phi0,'p=',p,'q=',q,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = lpi1_x_star(p,q,phi0,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = lpi1_epsilon_one(xstar,p,q,phi0)
           eps2 = lpi1_epsilon_two(xstar,p,q,phi0)
           eps3 = lpi1_epsilon_three(xstar,p,q,phi0)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('lpi1_predic.dat',p,q,phi0,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('lpi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/phi0,q,p/))

        end do

     end do

  enddo

  call aspicwrite_end()

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  p=4
  q=2
  phi0=0.1
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = lpi1_x_rrad(p,q,phi0,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = lpi1_epsilon_one(xstar,p,q,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = lpi1_x_endinf(p,q,phi0)
     eps1end =  lpi1_epsilon_one(xend,p,q,phi0)
     VendOverVstar = lpi1_norm_potential(xend,p,q,phi0)/lpi1_norm_potential(xstar,p,q,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = lpi1_x_rreh(p,q,phi0,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = lpi1_x_star(p,q,phi0,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program lpi1main
