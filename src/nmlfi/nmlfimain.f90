program nmlfimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar, HiggsCoupling
  use infinout, only : delete_file, livewrite

  use nmlficommon, only : nmlfi_norm_parametric_potential, nmlfi_norm_deriv_parametric_potential
  use nmlficommon, only : nmlfi_norm_deriv_second_parametric_potential
  use nmlficommon, only : nmlfi_parametric_epsilon_one, nmlfi_parametric_efold_primitive
  use nmlficommon, only : nmlfi_parametric_epsilon_two, nmlfi_parametric_epsilon_three
  use nmlficommon, only : nmlfi_deriv_x, nmlfi_deriv_second_x, nmlfi_parametric_hbar_endinf
  
  use nmlfisr, only : nmlfi_norm_potential, nmlfi_norm_deriv_potential, nmlfi_norm_deriv_second_potential
  use nmlfisr, only : nmlfi_epsilon_one, nmlfi_epsilon_two, nmlfi_epsilon_three, nmlfi_x_endinf  
  use nmlfisr, only : nmlfi_x_trajectory, nmlfi_x, nmlfi_hbar

  use srreheat, only : get_lnrreh_rrad, get_lnrreh_rhow, get_lnrrad_rhow
  use srreheat, only : ln_rho_reheat, ln_rho_endinf, log_energy_reheat_ingev
  use srreheat, only : potential_normalization
  use nmlfireheat, only : nmlfi_hbar_star, nmlfi_hbar_rrad, nmlfi_hbar_rreh, nmlfi_lnrhoreh_max
  use nmlfireheat, only : nmlfi_x_star, nmlfi_x_rrad, nmlfi_x_rreh

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  integer :: i,j,n,npts

  real(kp) :: f, fmin, fmax
  real(kp) :: lambda, ucte

  real(kp) :: hbar,hbarstar, hbarmin, hbarmax, hbarend
  real(kp) :: xend, xstar
  real(kp) :: xistar, xi
  
  real(kp) :: x, xmin, xmax
  real(kp) :: dx, d2x, d3x

  real(kp) :: V, dV, d2V, d3V
  real(kp) :: eps1, eps2, eps3

  real(kp) :: w, Pstar, ns, r, Treh, logErehGeV
  real(kp) :: lnRhoReh, lnRhoRehMin, lnRhoRehMax
  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, bfoldstar
  real(kp) :: M, Vstar, lnOmega4End
   

  lambda = nmlfiggsCoupling
  print *,'lambda=',lambda

  
  call delete_file('nmlfi_potential.dat')
  call delete_file('nmlfi_parametric_potential.dat')
  call delete_file('nmlfi_slowroll.dat')

  n=500

  xmin = 0._kp
  xmax = 10._kp

  xi = 20000._kp
  
  do i=1,n

     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)

     V = nmlfi_norm_potential(x,xi)

     call livewrite('nmlfi_potential.dat',x,V)

     eps1 = nmlfi_epsilon_one(x,xi)
     eps2 = nmlfi_epsilon_two(x,xi)
     eps3 = nmlfi_epsilon_three(x,xi)


     call livewrite('nmlfi_slowroll.dat',x,eps1,eps2,eps3)

  enddo


  do i=1,n

     hbar = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)
     x = nmlfi_x(hbar,xi)
     
     V = nmlfi_norm_parametric_potential(hbar,xi)

     call livewrite('nmlfi_parametric_potential.dat',hbar/sqrt(xi),V,x,hbar/sqrt(xi))

  
  enddo

  
  Pstar = powerAmpScalar

  call delete_file('nmlfi_prenmlfic.dat')
  call delete_file('nmlfi_nsr.dat')

  call aspicwrite_header('nmlfi',labeps12,labnsr,labbfoldreh)
 
  npts = 20

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = nmlfi_lnrhoreh_max(Pstar)

  print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     hbarstar = nmlfi_hbar_star(w,lnRhoReh,Pstar,bfoldstar,xistar)
    
     eps1 = nmlfi_parametric_epsilon_one(hbarstar,xistar)
     eps2 = nmlfi_parametric_epsilon_two(hbarstar,xistar)
     eps3 = nmlfi_parametric_epsilon_three(hbarstar,xistar)

     Vstar = nmlfi_norm_parametric_potential(hbarstar,xistar)
     M = potential_normalization(Pstar,eps1,Vstar)

     print *,'lnRhoReh= ',lnRhoReh, 'M= ', M, 'xistar/sqrt(lambda)= ',xistar/sqrt(lambda) &
          , 'bfoldstar= ',bfoldstar
     
     logErehGev = log_energy_reheat_ingev(lnRhoReh)
     Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

     ns = 1._kp-2._kp*eps1 - eps2
     r = 16._kp*eps1   

     call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/))

  enddo

  
  call aspicwrite_end()
  
! Test reheating with lnRrad and lnR

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'
  
  lnRradmin=-40
  lnRradmax = 10

  npts = 10

  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     hbarstar = nmlfi_hbar_rrad(lnRrad,Pstar,bfoldstar,xistar)

     print *,'lnRrad=',lnRrad, 'hbarstar=', hbarstar, 'xistar= ',xistar, 'bfoldstar= ',bfoldstar
     
     eps1 = nmlfi_parametric_epsilon_one(hbarstar,xistar)
     eps2 = nmlfi_parametric_epsilon_two(hbarstar,xistar)
     eps3 = nmlfi_parametric_epsilon_three(hbarstar,xistar)
     
!consistency test
!get lnR from lnRrad and check that it gives the same xstar
     hbarend = nmlfi_parametric_hbar_endinf(xistar)
     eps1end =  nmlfi_parametric_epsilon_one(hbarend,xistar)
     lnOmega4End = 2._kp*log(1._kp + hbarend*hbarend)

     VendOverVstar = nmlfi_norm_parametric_potential(hbarend,xistar) &
          /nmlfi_norm_parametric_potential(hbarstar,xistar)     
     
!in the Jordan Frame!!!     
     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar,lnOmega4End)
     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)

     hbarstar = nmlfi_hbar_rreh(lnR,Pstar)

     print *,'lnR',lnR,'hbarstar', hbarstar

!second consistency check
!get rhoreh for chosen w and check that hbarstar gotten tnmlfis way is the same
     w = 0._kp
!in the Jordan Frame!!!
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar,lnOmega4End)

     hbarstar = nmlfi_hbar_star(w,lnRhoReh,Pstar,xistar=xistar)
     xstar = nmlfi_x(hbarstar,xistar)

     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'hbarstar',hbarstar,'hbarend ',hbarend

     print *
     

    enddo



end program nmlfimain
