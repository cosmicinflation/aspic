program nmlfi1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar, HiggsCoupling
  use infinout, only : delete_file, livewrite

  use nmlficommon, only : nmlfi_parametric_epsilon_one, nmlfi_parametric_epsilon_two, nmlfi_parametric_epsilon_three
  use nmlficommon, only : nmlfi_x, nmlfi_norm_parametric_potential, pplus
  
  use nmlfi1sr, only : nmlfi1_norm_potential, nmlfi1_norm_deriv_potential, nmlfi1_norm_deriv_second_potential
  use nmlfi1sr, only : nmlfi1_epsilon_one, nmlfi1_epsilon_two, nmlfi1_epsilon_three, nmlfi1_x_endinf  
  use nmlfi1sr, only : nmlfi1_x_trajectory, nmlfi1_hbar_endinf, nmlfi1_ximax
  use nmlfi1sr, only : nmlfi1_numacc_efoldmax, nmlfi1_numacc_ximax

  use srreheat, only : get_lnrreh_rrad, get_lnrreh_rhow, get_lnrrad_rhow
  use srreheat, only : ln_rho_reheat, ln_rho_endinf, log_energy_reheat_ingev
  use srreheat, only : potential_normalization

  use nmlfi1reheat, only : nmlfi1_hbar_star, nmlfi1_hbar_rrad, nmlfi1_hbar_rreh, nmlfi1_parametric_lnrhoreh_max
  use nmlfi1reheat, only : nmlfi1_x_star, nmlfi1_x_rrad, nmlfi1_x_rreh

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  integer :: i,j,k,n,npts,nxi

  real(kp) :: p,V

  real(kp) :: hbar,hbarstar, hbarmin, hbarmax, hbarend
  real(kp) :: xend, xstar
  
  real(kp) :: x, xmin, xmax
  real(kp) :: eps1, eps2, eps3

  real(kp) :: w, Pstar, ns, r, Treh, logErehGeV
  real(kp) :: lnRhoReh, lnRhoRehMin, lnRhoRehMax
  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, bfoldstar
  real(kp) :: M, Vstar, lnOmega4End
   

  real(kp) :: xi, ximin,ximax, efoldMin, efoldMax
  
  call delete_file('nmlfi1_potential.dat')
  call delete_file('nmlfi1_parametric_potential.dat')
  call delete_file('nmlfi1_slowroll.dat')

  n=500

  xmin = 0.01_kp
  xmax = 10._kp

  xi = 10._kp
  p = 0.5
  
  do i=1,n

     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)

     V = nmlfi1_norm_potential(x,xi,p)

     call livewrite('nmlfi1_potential.dat',x,V)

     eps1 = nmlfi1_epsilon_one(x,xi,p)
     eps2 = nmlfi1_epsilon_two(x,xi,p)
     eps3 = nmlfi1_epsilon_three(x,xi,p)

     call livewrite('nmlfi1_slowroll.dat',x,eps1,eps2,eps3)

  enddo


  do i=1,n

     hbar = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)
     x = nmlfi_x(hbar,xi)
     
     V = nmlfi_norm_parametric_potential(hbar,xi,p)

     call livewrite('nmlfi1_parametric_potential.dat',hbar/sqrt(xi),V,x,hbar/sqrt(xi))

  
  enddo

  
  Pstar = powerAmpScalar
  w = 0._kp
  

  call aspicwrite_header('nmlfi1',labeps12,labnsr,labbfoldreh,(/'xi','p '/))
 
  npts = 20
  efoldMin = 120._kp
  
  lnRhoRehMin = lnRhoNuc

  do k=1,8

     p = k

     nxi=100
     
     ximin = 1d-6
     ximax = 1d6

     if (p.ge.pplus) then
        print *,'ximaxth= ',nmlfi1_ximax(p)
        ximax = min(ximax,nmlfi1_ximax(p))
     endif
     
     ximax = min(ximax,nmlfi1_numacc_ximax(efoldMin,p))
     print *,'ximin, ximax= ',ximin,ximax
     
     do j=1,nxi
        
        xi = exp(log(ximin) + real(j-1,kp)*(log(ximax)-log(ximin))/real(nxi-1,kp))

        efoldMax = nmlfi1_numacc_efoldmax(xi,p)

        if (efoldMax.lt.efoldMin) then
           print *,'not enough e-folds: xi= efoldMax= ',xi,efoldMax
           cycle
        end if
        
        hbarend = nmlfi1_hbar_endinf(xi,p)

        lnRhoRehMax = nmlfi1_parametric_lnrhoreh_max(xi,p,hbarend,Pstar)

        print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           hbarstar = nmlfi1_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)

!           print *,'test',hbarstar,hbarend,xi,p,w

           
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

  xi = 0.1
  p=5
  npts = 10
  
  hbarend = nmlfi1_hbar_endinf(xi,p)
  
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     hbarstar = nmlfi1_hbar_rrad(xi,p,hbarend,lnRrad,Pstar,bfoldstar)

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

     hbarstar = nmlfi1_hbar_rreh(xi,p,hbarend,lnR)

     print *,'lnR',lnR,'hbarstar', hbarstar

!second consistency check
!get rhoreh for chosen w and check that hbarstar gotten tnmlfi1s way is the same
     w = 0._kp
!in the Jordan Frame!!!
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar,lnOmega4End)

     hbarstar = nmlfi1_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar)
     xstar = nmlfi_x(hbarstar,xi)
     
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'hbarstar',hbarstar,'hbarend ',hbarend

     print *
     

    enddo



end program nmlfi1main
