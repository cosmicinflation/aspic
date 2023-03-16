program nmlfi3main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar, HiggsCoupling
  use infinout, only : delete_file, livewrite

  use nmlficommon, only : pplus, pminus, nmlfi_xizero
  use nmlficommon, only : nmlfi_parametric_epsilon_one, nmlfi_parametric_epsilon_two, nmlfi_parametric_epsilon_three
  use nmlficommon, only : nmlfi_parametric_potential, nmlfi_x
  
  use nmlfi3sr, only : nmlfi3_norm_potential, nmlfi3_norm_deriv_potential, nmlfi3_norm_deriv_second_potential
  use nmlfi3sr, only : nmlfi3_epsilon_one, nmlfi3_epsilon_two, nmlfi3_epsilon_three, nmlfi3_x_endinf  
  use nmlfi3sr, only : nmlfi3_x_trajectory

  use srreheat, only : get_lnrreh_rrad, get_lnrreh_rhow, get_lnrrad_rhow
  use srreheat, only : ln_rho_reheat, ln_rho_endinf, log_energy_reheat_ingev
  use srreheat, only : potential_normalization
  use nmlfi3reheat, only : nmlfi3_hbar_star, nmlfi3_hbar_rrad, nmlfi3_hbar_rreh, nmlfi3_lnrhoreh_max
  use nmlfi3reheat, only : nmlfi3_x_star, nmlfi3_x_rrad, nmlfi3_x_rreh

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  integer :: i,j,k,n,npts,nxi,nend

  real(kp) :: hbar,hbarstar, hbarmin, hbarmax, hbarend
  real(kp) :: xstar
  
  real(kp) :: x, xmin, xmax

  real(kp) :: eps1, eps2, eps3

  real(kp) :: w, Pstar, ns, r, Treh, logErehGeV
  real(kp) :: lnRhoReh, lnRhoRehMin, lnRhoRehMax
  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, bfoldstar
  real(kp) :: M, Vstar, lnOmega4End
   

  real(kp) :: xi, ximin,ximax
  real(kp) :: p, pmin,pmax
  real(kp) :: xend, xendmin, xendmax
  
  
  Pstar = powerAmpScalar
  

  call aspicwrite_header('nmlfi3',labeps12,labnsr,labbfoldreh,(/'xi','p '/))
 
  npts = 20

  np = 5
  pmin = 0.1
  pmax = pminus
  
  nxi=100
  
  lnRhoRehMin = lnRhoNuc

  do k=1,np

     p = pmin +real(k-1,kp)*(pmax-pmin)/real(np-1,kp)
     
     do j=1,nxi

        ximin = nmlfi_xizero(p)
        ximax = 100._kp * xmin
        
        xi = ximin + real(j-1,kp)*(ximax-ximin)/real(nxi-1,kp)
    
        hbarend = nmlfi3_hbar_endinf(xi,p)
        lnRhoRehMax = nmlfi3_lnrhoreh_max(xi,p,hbarend,Pstar)

        print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           hbarstar = nmlfi3_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)

           eps1 = nmlfi_parametric_epsilon_one(hbarstar,xi,p)
           eps2 = nmlfi_parametric_epsilon_two(hbarstar,xi,p)
           eps3 = nmlfi_parametric_epsilon_three(hbarstar,xi,p)

           Vstar = nmlfi_norm_parametric_potential(hbarstar,xi,p)
           M = potential_normalization(Pstar,eps1,Vstar)

           print *,'lnRhoReh= ',lnRhoReh, 'M= ', M, 'xi= ',xi, 'p= ',p &
                , 'bfoldstar= ',bfoldstar

           logErehGev = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp-2._kp*eps1 - eps2
           r = 16._kp*eps1   

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xi,p/))

        enddo

     enddo
  enddo
  
  call aspicwrite_end()
  
! Test reheating with lnRrad and lnR

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'
  
  lnRradmin=-40
  lnRradmax = 10

  p = 0.5_kp*pminus
  xi = 10._kp * nmlfi_xizero(p)
  npts = 10
  
  hbarend = nmlfi3_parametric_hbar_endinf(xi,p)
  
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     hbarstar = nmlfi3_hbar_rrad(xi,p,hbarend,lnRrad,Pstar,bfoldstar)

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

     hbarstar = nmlfi3_hbar_rreh(xi,p,hbarend,lnR)

     print *,'lnR',lnR,'hbarstar', hbarstar

!second consistency check
!get rhoreh for chosen w and check that hbarstar gotten tnmlfi3s way is the same
     w = 0._kp
!in the Jordan Frame!!!
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar,lnOmega4End)

     hbarstar = nmlfi3_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar)
     xstar = nmlfi_x(hbarstar,xi)

     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'hbarstar',hbarstar,'hbarend ',hbarend

     print *
     

    enddo



end program nmlfi3main
