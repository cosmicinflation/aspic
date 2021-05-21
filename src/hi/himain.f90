program himain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use infinout, only : delete_file, livewrite

  use hicommon, only : hi_norm_parametric_potential, hi_norm_deriv_parametric_potential
  use hicommon, only : hi_norm_deriv_second_parametric_potential
  use hicommon, only : hi_parametric_epsilon_one, hi_parametric_efold_primitive
  use hicommon, only : hi_parametric_epsilon_two, hi_parametric_epsilon_three
  use hicommon, only : hi_deriv_x, hi_deriv_second_x, hi_deriv_third_x 
  
  use hisr, only : hi_norm_potential, hi_norm_deriv_potential, hi_norm_deriv_second_potential
  use hisr, only : hi_epsilon_one, hi_epsilon_two, hi_epsilon_three, hi_x_enhinf  
  use hisr, only : hi_x_trajectory

  use srreheat, only : get_lnrreh_rrad, get_lnrreh_rhow, get_lnrrad_rhow
  use srreheat, only : ln_rho_reheat, ln_rho_enhinf, log_energy_reheat_ingev
  use hireheat, only : hi_hbar_star, hi_hbar_rrad, hi_hbar_rreh, hi_lnrhoreh_max
  use hireheat, only : hi_x_star, hi_x_rrad, hi_x_rreh

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  integer :: i,j,n,npts

  real(kp) :: f, fmin, fmax
  real(kp) :: lambda, ucte

  real(kp) :: hbar, hbarmin, hbarmax, hbarend
  real(kp) :: xend, xstar
  real(kp) :: xistar
  
  real(kp) :: x, xmin, xmax
  real(kp) :: dx, d2x, d3x

  real(kp) :: V, dV, d2V, d3V
  real(kp) :: eps1, eps2, eps3

  real(kp) :: w, Pstar, ns, r, Treh, logErehGeV
  real(kp) :: lnRhoReh, lnRhoRehMin, lnRhoRehMax
  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, bfoldstar
   
  

  call delete_file('hi_potential.dat')
  call delete_file('hi_slowroll.dat')

  n=250

  xmin = 0.0_kp
  xmax = 10._kp

  do i=1,n

     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)

     V = hi_norm_potential(x)

     call livewrite('hi_potential.dat',x,V)

     eps1 = hi_epsilon_one(x)
     eps2 = hi_epsilon_two(x)
     eps3 = hi_epsilon_three(x)


     call livewrite('si_slowroll.dat',x,eps1,eps2,eps3)

  enddo

  Pstar = powerAmpScalar



!  f=1e-3
!  print *,hi_lnrhoreh_max(f,Pstar)
!  stop

  call delete_file('hi_prehic.dat')
  call delete_file('hi_nsr.dat')

  call aspicwrite_header('hi',labeps12,labnsr,labbfoldreh)
 
  npts = 20

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = hi_lnrhoreh_max(Pstar)

  print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     hbarstar = hi_hbar_star(w,lnRhoReh,Pstar,bfoldstar,xistar)

     print *,'lnRhoReh= ',lnRhoReh, 'xistar= ',xistar, 'bfoldstar= ',bfoldstar

     eps1 = hi_parametric_epsilon_one(hbarstar,xistar)
     eps2 = hi_parametric_epsilon_two(hbarstar,xistar)
     eps3 = hi_parametric_epsilon_three(hbarstar,xistar)

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

     hbarstar = hi_hbar_rrad(f,lnRrad,Pstar,bfoldstar,xistar)



     print *,'lnRrad=',lnRrad, 'hbarstar=', hbarstar, 'xistar= ',xistar, 'bfoldstar= ',bfoldstar
     
     eps1 = hi_parametric_epsilon_one(hbarstar,xistar)
     eps2 = hi_parametric_epsilon_two(hbarstar,xistar)
     eps3 = hi_parametric_epsilon_three(hbarstar,xistar)
     
!consistency test
!get lnR from lnRrad and check that it gives the same xstar
     hbarend = hi_parametric_hbar_endinf(xistar)
     eps1end =  hi_parametric_epsilon_one(hbarend,xistar)
     VendOverVstar = hi_norm_parametric_potential(hbarend,xistar) &
          /hi_norm_parametric_potential(hbarstar,xistar)

     lnRhoEnd = ln_rho_enhinf(Pstar,eps1,eps1End,VendOverVstar)
     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)

     hbarstar = hi_hbar_rreh(lnR,Pstar)

     print *,'lnR',lnR,'hbarstar', hbarstar

!second consistency check
!get rhoreh for chosen w and check that hbarstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     hbarstar = hi_hbar_star(w,lnRhoReh,Pstar)
     xstar = hi_x(hbarstar)

     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'hbarstar',hbarstar



    enddo



end program himain
