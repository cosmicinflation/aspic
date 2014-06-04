program nfi1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use nfi1sr, only : nfi1_epsilon_one, nfi1_epsilon_two, nfi1_epsilon_three
  use nfi1reheat, only : nfi1_lnrhoreh_max, nfi1_lnrhoreh_fromepsilon 
  use nfi1reheat, only : nfi1_xp_fromepsilon, nfi1_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  use srflow, only : scalar_spectral_index, tensor_to_scalar_ratio

  use nfi1sr, only : nfi1_norm_potential, nfi1_x_endinf
  use nfi1reheat, only : nfi1_x_rreh, nfi1_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use srreheat, only : potential_normalization, primscalar

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i
  integer :: npts = 20

  real(kp) :: a,b,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: M, As

  Pstar = powerAmpScalar

  call delete_file('nfi1_predic.dat')
  call delete_file('nfi1_nsr.dat')

  a = 0._kp 
  b = 1.1_kp

  w = 0
 
 
    lnRhoRehMin = lnRhoNuc
    lnRhoRehMax = nfi1_lnrhoreh_max(p,Pstar)

    print *,'a= b= ',a,b,'lnRhoRehMin=',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

    do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = nfi1_x_star(a,b,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar

       eps1 = nfi1_epsilon_one(xstar,a,b)
       eps2 = nfi1_epsilon_two(xstar,a,b)
       eps3 = nfi1_epsilon_three(xstar,a,b)       

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('nfi1_predic.dat',a,b,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('nfi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

  enddo

! Test reheating with lnRrad and lnR

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'
  
  lnRradmin=-42
  lnRradmax = 10
  a = 0.1
  b = 1.3

  do i=1,npts

       lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = nfi1_x_rrad(p,lnRrad,Pstar,bfoldstar)

       print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

       eps1 = nfi1_epsilon_one(xstar,a,b)
       eps2 = nfi1_epsilon_two(xstar,a,b)
       eps3 = nfi1_epsilon_three(xstar,a,b)
       
!consistency test
!get lnR from lnRrad and check that it gives the same xstar
       xend = nfi1_x_endinf(a,b)
       eps1end =  nfi1_epsilon_one(xend,a,b)
       VendOverVstar = nfi1_norm_potential(xend,a,b) &
            /nfi1_norm_potential(xstar,a,b)

       lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

       lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
       xstar = nfi1_x_rreh(a,b,lnR,bfoldstar)
       print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar


!second consistency check
!get rhoreh for chosen w and check that xstar gotten this way is the same
       w = 0._kp
       lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)
       
       xstar = nfi1_x_star(a,b,w,lnRhoReh,Pstar,bfoldstar)
       print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
            ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

!check potential normalization/Pstar
       M = potential_normalization(Pstar,(/eps1,eps2,eps3/) &
            ,nfi1_norm_potential(xstar,a,b))
       As = primscalar(M,(/eps1,eps2,eps3/),nfi1_norm_potential(xstar,a,b))
       print *,'M=',M,As,Pstar
                    
    enddo


  end program nfi1main
