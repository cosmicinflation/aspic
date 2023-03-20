program nmlfi2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar, HiggsCoupling
  use infinout, only : delete_file, livewrite

  use nmlficommon, only : pplus, pminus, nmlfi_xizero
  use nmlficommon, only : nmlfi_parametric_epsilon_one, nmlfi_parametric_epsilon_two, nmlfi_parametric_epsilon_three
  use nmlficommon, only : nmlfi_norm_parametric_potential, nmlfi_x, nmlfi_hbar_potmax
  
  use nmlfi2sr, only : nmlfi2_norm_potential, nmlfi2_norm_deriv_potential, nmlfi2_norm_deriv_second_potential
  use nmlfi2sr, only : nmlfi2_epsilon_one, nmlfi2_epsilon_two, nmlfi2_epsilon_three, nmlfi2_x_endinf  
  use nmlfi2sr, only : nmlfi2_x_trajectory, nmlfi2_hbar_endinf
  use nmlfi2sr, only : nmlfi2_numacc_efoldmax, nmlfi2_numacc_hbarinimin

  use srreheat, only : get_lnrreh_rrad, get_lnrreh_rhow, get_lnrrad_rhow
  use srreheat, only : ln_rho_reheat, ln_rho_endinf, log_energy_reheat_ingev
  use srreheat, only : potential_normalization
  use nmlfi2reheat, only : nmlfi2_hbar_star, nmlfi2_hbar_rrad, nmlfi2_hbar_rreh
  use nmlfi2reheat, only : nmlfi2_x_star, nmlfi2_x_rrad, nmlfi2_x_rreh, nmlfi2_parametric_lnrhoreh_max

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  integer :: i,j,k,n,np,npts,nxi

  real(kp) :: hbar,hbarstar, hbarmin, hbarmax, hbarend
  real(kp) :: xend, xstar
  
  real(kp) :: x, xmin, xmax, efoldMax

  real(kp) :: eps1, eps2, eps3

  real(kp) :: w, Pstar, ns, r, Treh, logErehGeV
  real(kp) :: lnRhoReh, lnRhoRehMin, lnRhoRehMax
  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, bfoldstar
  real(kp) :: M, Vstar, lnOmega4End
   

  real(kp) :: xi, ximin,ximax
  real(kp) :: p, pmin,pmax, xizero
  
  
  Pstar = powerAmpScalar
  w = 0._kp

  call aspicwrite_header('nmlfi2',labeps12,labnsr,labbfoldreh,(/'xi','p '/))
 
  npts = 20

  np = 5
  pmin = 0.01_kp
  pmax = 0.1_kp
  
  nxi=10
  
  lnRhoRehMin = lnRhoNuc

  do k=1,np

     p = pmin +real(k-1,kp)*(pmax-pmin)/real(np-1,kp)

     xizero = nmlfi_xizero(p)

     print *,'p= xizero= ',p,xizero
     
     do j=1,nxi
        
        ximin = 1.01_kp * xizero
        ximax = 1d3
        
        xi = exp(log(ximin) + real(j-1,kp)*(log(ximax)-log(ximin))/real(nxi-1,kp))

        efoldMax = nmlfi2_numacc_efoldmax(xi,p)
        
        print *, 'efoldMax=',efoldMax

        if (efoldMax.lt.120._kp) then
           write(*,*)'not enough e-folds to be viable!'
           cycle
        endif
        hbarend = nmlfi2_hbar_endinf(xi,p)

        lnRhoRehMax = nmlfi2_parametric_lnrhoreh_max(xi,p,hbarend,Pstar)

        print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax
        print *, 'efoldmax= ',efoldMAx

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           hbarstar = nmlfi2_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)

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
  
  lnRradmin=-46
  lnRradmax = -40

  p = 0.1*pminus
  xi = 1.1_kp * nmlfi_xizero(p)
  npts = 10
  
  hbarend = nmlfi2_hbar_endinf(xi,p)
  
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     hbarstar = nmlfi2_hbar_rrad(xi,p,hbarend,lnRrad,Pstar,bfoldstar)

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

     hbarstar = nmlfi2_hbar_rreh(xi,p,hbarend,lnR)

     print *,'lnR',lnR,'hbarstar', hbarstar

!second consistency check
!get rhoreh for chosen w and check that hbarstar gotten tnmlfi2s way is the same
     w = 0._kp
!in the Jordan Frame!!!
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar,lnOmega4End)

     hbarstar = nmlfi2_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar)
     xstar = nmlfi_x(hbarstar,xi)

     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'hbarstar',hbarstar,'hbarend ',hbarend

     print *
     

    enddo



end program nmlfi2main
