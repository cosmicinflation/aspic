!test the reheating derivation from slow-roll
program psnimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use psnisr, only : psni_epsilon_one, psni_epsilon_two, psni_epsilon_three
  use psnireheat, only : psni_lnrhoreh_max, psni_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use psnisr, only : psni_norm_potential, psni_x_endinf
  use psnireheat, only : psni_x_rreh, psni_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nf

  real(kp) :: alpha,f,w,bfoldstar,alphamin,alphamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r
  real(kp), dimension(:), allocatable :: fvalues


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
  call delete_file('psni_predic.dat')
  call delete_file('psni_nsr.dat')

  npts = 20 
  nalpha=200
  w=0._kp
  !  w = 1._kp/3._kp


  nf=3
  allocate(fvalues(1:3))
  fvalues(3)=0.001_kp
  fvalues(2)=0.1_kp
  fvalues(1)=10._kp

  do j=1,nf
     f=fvalues(j)

     alphamax=10._kp**(-1._kp)*f**2
     alphamin=10._kp**(-8._kp)*f**2

     do k=1,nalpha
        alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !log step

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = psni_lnrhoreh_max(alpha,f,Pstar)

        print *,'alpha=',alpha,'f=',f,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = psni_x_star(alpha,f,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = psni_epsilon_one(xstar,alpha,f)
           eps2 = psni_epsilon_two(xstar,alpha,f)
           eps3 = psni_epsilon_three(xstar,alpha,f)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('psni_predic.dat',alpha,f,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('psni_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 1e-4
  f=0.1
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = psni_x_rrad(alpha,f,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = psni_epsilon_one(xstar,alpha,f)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = psni_x_endinf(alpha,f)
     eps1end =  psni_epsilon_one(xend,alpha,f)
     VendOverVstar = psni_norm_potential(xend,alpha,f)/psni_norm_potential(xstar,alpha,f)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = psni_x_rreh(alpha,f,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = psni_x_star(alpha,f,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program psnimain
