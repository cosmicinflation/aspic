program nmlfi3main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar, HiggsCoupling
  use infinout, only : delete_file, livewrite

  use nmlficommon, only : pplus, pminus, nmlfi_xizero
  use nmlficommon, only : nmlfi_parametric_epsilon_one, nmlfi_parametric_epsilon_two, nmlfi_parametric_epsilon_three
  use nmlficommon, only : nmlfi_norm_parametric_potential, nmlfi_x, nmlfi_hbar_potmax
  use nmlficomreh, only : nmlfi_gravity_mass_scale
  
  use nmlfi3sr, only : nmlfi3_norm_potential, nmlfi3_norm_deriv_potential, nmlfi3_norm_deriv_second_potential
  use nmlfi3sr, only : nmlfi3_epsilon_one, nmlfi3_epsilon_two, nmlfi3_epsilon_three
  use nmlfi3sr, only : nmlfi3_x_trajectory, nmlfi3_numacc_hbarinimin, nmlfi3_numacc_efoldmax
  use nmlfi3sr, only : nmlfi3_numacc_hbarendmin
  
  use srreheat, only : get_lnrreh_rrad, get_lnrreh_rhow, get_lnrrad_rhow
  use srreheat, only : ln_rho_reheat, ln_rho_endinf, log_energy_reheat_ingev
  use srreheat, only : potential_normalization
  use nmlfi3reheat, only : nmlfi3_hbar_star, nmlfi3_hbar_rrad, nmlfi3_hbar_rreh
  use nmlfi3reheat, only : nmlfi3_x_star, nmlfi3_x_rrad, nmlfi3_x_rreh,nmlfi3_parametric_lnrhoreh_max

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  integer :: i,j,k,l,n,npts,nend

  real(kp) :: hbar,hbarstar, hbarmin, hbarmax
  real(kp) :: hbarend, hbarendmax, hbarendmin
  real(kp) :: xstar
  
  real(kp) :: x, xmin, xmax

  real(kp) :: eps1, eps2, eps3

  real(kp) :: w, Pstar, ns, r, Treh, logErehGeV
  real(kp) :: lnRhoReh, lnRhoRehMin, lnRhoRehMax
  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, bfoldstar
  real(kp) :: M, Mg, Vstar, lnOmega4End
   

  real(kp) :: xi, ximin,ximax,xizero
  real(kp) :: p, pmin,pmax, efoldMin
  real(kp) :: xend, xendmin, xendmax

  integer, parameter :: np = 3
  integer, parameter :: nxi = 3
  real(kp), dimension(nxi,np) :: xivec
  
  
  Pstar = powerAmpScalar
  w = 0._kp
  efoldMin = 120._kp
  
  call aspicwrite_header('nmlfi3s',labeps12,labnsr,labbfoldreh,(/'chiend','xi    ','p     '/))
 
  npts = 20

  pmin = 0.1_kp
  pmax = 0.4_kp

  
  xivec(1,1) = 0.001_kp
  xivec(2,1) = 0.0025_kp
  xivec(3,1) = 0.004_kp

  xivec(1,2) = 0.001_kp
  xivec(2,2) = 0.003_kp
  xivec(3,2) = 0.004_kp

  xivec(1,3) = 0.001_kp
  xivec(2,3) = 0.003_kp
  xivec(3,3) = 0.004_kp
  
!reminder  
  if (pmax.gt.pminus) stop 'nmlfi3s is for p<p- and xi<xizero'
  

  nend = 1000

  
  lnRhoRehMin = lnRhoNuc

  do k=1,np

     p = pmin +real(k-1,kp)*(pmax-pmin)/real(np-1,kp)
     xizero = nmlfi_xizero(p)
     
     !nmlfi3 with p<p- and xi<xizero, eps1 is always < 1. If slow-roll has to be enforced, xi cannot be too close of xizero

     do j=1,nxi

        xi = xivec(j,k)

        hbarendmin = nmlfi3_numacc_hbarendmin(efoldMin,xi,p)
        hbarendmax = 10._kp * hbarendmin
        print *
        print *,'hbarendmin= hbarendmax= ',hbarendmin,hbarendmax
        print *,'p= pminus= ',p,pminus
        print *,'xi= xizero(p)= ',xi,xizero

        if (xi.gt.xizero) stop 'buffoon!'

        do l=1,nend

           hbarend = hbarendmin + real(l-1,kp)*(hbarendmax-hbarendmin)/real(nend-1,kp)

           xend = nmlfi_x(hbarend,xi)

           lnRhoRehMax = nmlfi3_parametric_lnrhoreh_max(xi,p,hbarend,Pstar)

           print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

           do i=1,npts

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              hbarstar = nmlfi3_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)

              eps1 = nmlfi_parametric_epsilon_one(hbarstar,xi,p)
              eps2 = nmlfi_parametric_epsilon_two(hbarstar,xi,p)
              eps3 = nmlfi_parametric_epsilon_three(hbarstar,xi,p)

              Vstar = nmlfi_norm_parametric_potential(hbarstar,xi,p)
              M = potential_normalization(Pstar,eps1,Vstar)
              Mg = nmlfi_gravity_mass_scale(xi,p,hbarend)

              print *,'M= Mg= ',M,Mg

              print *,'lnRhoReh= ',lnRhoReh, 'xi= ',xi, 'p= ',p &
                   , 'bfoldstar= ',bfoldstar

              logErehGev = log_energy_reheat_ingev(lnRhoReh)
              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

              ns = 1._kp-2._kp*eps1 - eps2
              r = 16._kp*eps1   

              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend,xi,p/))

           enddo

        enddo

     enddo

  enddo
  
  call aspicwrite_end()

  call aspicwrite_header('nmlfi3l',labeps12,labnsr,labbfoldreh,(/'chiend','xi    ','p     '/))
 
  npts = 20

  pmin = 0.6_kp
  pmax = 3.5_kp

  if ((pmin.le.pminus).or.(p.ge.4._kp)) stop 'nmlfi3l is for p>pminus'


  xivec(1,1) = 0.001_kp
  xivec(2,1) = 0.003_kp
  xivec(3,1) = 0.004_kp

  xivec(1,2) = 0.002_kp
  xivec(2,2) = 0.004_kp
  xivec(3,2) = 0.005_kp

  xivec(1,3) = 0.01_kp
  xivec(2,3) = 0.02_kp
  xivec(3,3) = 0.03_kp
  
  
  nend = 1000
  
  lnRhoRehMin = lnRhoNuc

  do k=1,np

     p = pmin +real(k-1,kp)*(pmax-pmin)/real(np-1,kp)
     
     do j=1,nxi
        
        xi = xivec(j,k)


        hbarendmin = nmlfi3_numacc_hbarendmin(efoldMin,xi,p)
        hbarendmax = 10._kp * hbarendmin

        print *
        print *,'hbarendmin= hbarendmax= ',hbarendmin,hbarendmax
        print *,'p= pminus= ',p,pminus
        print *,'xi= ',xi
        
        
        do l=1,nend

           hbarend = hbarendmin + real(l-1,kp)*(hbarendmax-hbarendmin)/real(nend-1,kp)

           xend = nmlfi_x(hbarend,xi)

           lnRhoRehMax = nmlfi3_parametric_lnrhoreh_max(xi,p,hbarend,Pstar)

           print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

           do i=1,npts

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              hbarstar = nmlfi3_hbar_star(xi,p,hbarend,w,lnRhoReh,Pstar,bfoldstar)

              eps1 = nmlfi_parametric_epsilon_one(hbarstar,xi,p)
              eps2 = nmlfi_parametric_epsilon_two(hbarstar,xi,p)
              eps3 = nmlfi_parametric_epsilon_three(hbarstar,xi,p)

              Vstar = nmlfi_norm_parametric_potential(hbarstar,xi,p)
              M = potential_normalization(Pstar,eps1,Vstar)
              Mg = nmlfi_gravity_mass_scale(xi,p,hbarend)

              print *,'M= Mg= ',M,Mg
              
              print *,'lnRhoReh= ',lnRhoReh,'xi= ',xi, 'p= ',p &
                   , 'bfoldstar= ',bfoldstar

              logErehGev = log_energy_reheat_ingev(lnRhoReh)
              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

              ns = 1._kp-2._kp*eps1 - eps2
              r = 16._kp*eps1   

              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend,xi,p/))

           enddo

        enddo

     enddo

  enddo
  
  call aspicwrite_end()

  
! Test reheating with lnRrad and lnR

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'
  
  lnRradmin=-40
  lnRradmax = 10

  p = 3._kp
  xi = 10._kp
  npts = 10
  
  hbarend = 1000._kp*nmlfi_hbar_potmax(xi,p)
  
  
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
