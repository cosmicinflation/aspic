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


  integer :: i,j,k,n,npts

  real(kp) :: p

  real(kp) :: hbar,hbarstar, hbarmin, hbarmax, hbarend
  real(kp) :: xend, xstar
  real(kp) :: xi
  
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
   

  real(kp) :: xi, ximin,ximax
  
  call delete_file('nmlfi_potential.dat')
  call delete_file('nmlfi_parametric_potential.dat')
  call delete_file('nmlfi_slowroll.dat')

  n=500

  xmin = 0._kp
  xmax = 10._kp

  xi = 20000._kp
  p = 2
  
  do i=1,n

     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)

     V = nmlfi_norm_potential(x,xi,p)

     call livewrite('nmlfi_potential.dat',x,V)

     eps1 = nmlfi_epsilon_one(x,xi,p)
     eps2 = nmlfi_epsilon_two(x,xi,p)
     eps3 = nmlfi_epsilon_three(x,xi,p)

     call livewrite('nmlfi_slowroll.dat',x,eps1,eps2,eps3)

  enddo


  do i=1,n

     hbar = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)
     x = nmlfi_x(hbar,xi)
     
     V = nmlfi_norm_parametric_potential(hbar,xi,p)

     call livewrite('nmlfi_parametric_potential.dat',hbar/sqrt(xi),V,x,hbar/sqrt(xi))

  
  enddo

  
  Pstar = powerAmpScalar
  

  call aspicwrite_header('nmlfi',labeps12,labnsr,labbfoldreh,(/'xi','p '/)))
 
  npts = 20
  ximin = 1d-3
  ximax = 1d3
  
  lnRhoRehMin = lnRhoNuc

  do k=2,8

     p = k
     
     do j=1,nxi
        xi = ximin + real(j-1,kp)*(ximax-ximin)/real(j-1,kp)
  
  
        hbarend = nmlfi_parametric_hbar_endinf(xi,p)
        lnRhoRehMax = nmlfi_lnrhoreh_max(xi,p,hbarend,Pstar)

        print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           hbarstar = nmlfi_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)

           eps1 = nmlfi_parametric_epsilon_one(hbarstar,xi,p)
           eps2 = nmlfi_parametric_epsilon_two(hbarstar,xi,p)
           eps3 = nmlfi_parametric_epsilon_three(hbarstar,xi,p)

           Vstar = nmlfi_norm_parametric_potential(hbarstar,xi,p)
           M = potential_normalization(Pstar,eps1,Vstar)

           print *,'lnRhoReh= ',lnRhoReh, 'M= ', M, 'xi= ',xistar, 'p= ',p &
                , 'bfoldstar= ',bfoldstar

           logErehGev = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp-2._kp*eps1 - eps2
           r = 16._kp*eps1   

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/))

        enddo

     enddo
  enddo
  
  call aspicwrite_end()
  
! Test reheating with lnRrad and lnR

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'
  
  lnRradmin=-40
  lnRradmax = 10

  xi = 1000
  p=2
  npts = 10
  
  hbarend = nmlfi_parametric_hbar_endinf(xi,p)
  
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     hbarstar = nmlfi_hbar_rrad(xi,p,hbarend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad, 'hbarstar=', hbarstar, 'xi= ',xi, 'p= ',p, 'bfoldstar= ',bfoldstar
     
     eps1 = nmlfi_parametric_epsilon_one(hbarstar,xi,p)
     eps2 = nmlfi_parametric_epsilon_two(hbarstar,xi,p)
     eps3 = nmlfi_parametric_epsilon_three(hbarstar,xi,p)
     
!consistency test
!get lnR from lnRrad and check that it gives the same xstar
     eps1end =  nmlfi_parametric_epsilon_one(hbarend,xi,p)
     lnOmega4End = 2._kp*log(1._kp + hbarend*hbarend)

     VendOverVstar = nmlfi_norm_parametric_potential(hbarend,xi,p) &
          /nmlfi_norm_parametric_potential(hbarstar,xi,p)     
     
!in the Jordan Frame!!!     
     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar,lnOmega4End)
     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)

     hbarstar = nmlfi_hbar_rreh(xi,p,hbarend,lnRreh)

     print *,'lnR',lnR,'hbarstar', hbarstar

!second consistency check
!get rhoreh for chosen w and check that hbarstar gotten tnmlfis way is the same
     w = 0._kp
!in the Jordan Frame!!!
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar,lnOmega4End)

     hbarstar = nmlfi_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar)
     xstar = nmlfi_x(hbarstar,xi)

     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'hbarstar',hbarstar,'hbarend ',hbarend

     print *
     

    enddo



end program nmlfimain
