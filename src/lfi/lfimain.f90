!test the reheating derivation from slow-roll
program lfimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lfisr, only : lfi_epsilon_one, lfi_epsilon_two, lfi_epsilon_three
  use lfireheat, only : lfi_lnrhoreh_max, lfi_lnrhoreh_fromepsilon 
  use lfireheat, only : lfi_xp_fromepsilon, lfi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  use srflow, only : scalar_spectral_index, tensor_to_scalar_ratio

  use lfisr, only : lfi_norm_potential, lfi_x_endinf
  use lfireheat, only : lfi_x_rreh, lfi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i
  integer :: npts = 20

  real(kp) :: p,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar

  call delete_file('lfi_predic.dat')
  call delete_file('lfi_nsr.dat')

  p = 0_kp 
  do while (p<13._kp)
    
     p=p+1_kp
     w = (p-2)/(p+2)
 
    lnRhoRehMin = lnRhoNuc
    lnRhoRehMax = lfi_lnrhoreh_max(p,Pstar)

    print *,'p=',p,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

    do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = lfi_x_star(p,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar

       eps1 = lfi_epsilon_one(xstar,p)
       eps2 = lfi_epsilon_two(xstar,p)
       eps3 = lfi_epsilon_three(xstar,p)

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('lfi_predic.dat',p,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('lfi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('lfi2_predic_summarized.dat')
  p=2._kp
  w = (p-2)/(p+2)
  lnRhoReh = lnRhoNuc
  xstar = lfi_x_star(p,w,lnRhoReh,Pstar,bfoldstar)
  eps1 = lfi_epsilon_one(xstar,p)
  eps2 = lfi_epsilon_two(xstar,p)
  eps3 = lfi_epsilon_three(xstar,p)
  ns = 1._kp - 2._kp*eps1 - eps2
  r =16._kp*eps1
  call livewrite('lfi2_predic_summarized.dat',eps1,eps2,eps3,r,ns)
  lnRhoReh = lfi_lnrhoreh_max(p,Pstar)
  xstar = lfi_x_star(p,w,lnRhoReh,Pstar,bfoldstar)
  eps1 = lfi_epsilon_one(xstar,p)
  eps2 = lfi_epsilon_two(xstar,p)
  eps3 = lfi_epsilon_three(xstar,p)
  ns = 1._kp - 2._kp*eps1 - eps2
  r =16._kp*eps1
  call livewrite('lfi2_predic_summarized.dat',eps1,eps2,eps3,r,ns)

  call delete_file('lfi4_predic_summarized.dat')
  p=4._kp
  w = (p-2)/(p+2)
  lnRhoReh = lnRhoNuc
  xstar = lfi_x_star(p,w,lnRhoReh,Pstar,bfoldstar)
  eps1 = lfi_epsilon_one(xstar,p)
  eps2 = lfi_epsilon_two(xstar,p)
  eps3 = lfi_epsilon_three(xstar,p)
  ns = 1._kp - 2._kp*eps1 - eps2
  r =16._kp*eps1
  call livewrite('lfi4_predic_summarized.dat',eps1,eps2,eps3,r,ns)
  lnRhoReh = lfi_lnrhoreh_max(p,Pstar)
  xstar = lfi_x_star(p,w,lnRhoReh,Pstar,bfoldstar)
  eps1 = lfi_epsilon_one(xstar,p)
  eps2 = lfi_epsilon_two(xstar,p)
  eps3 = lfi_epsilon_three(xstar,p)
  ns = 1._kp - 2._kp*eps1 - eps2
  r =16._kp*eps1
  call livewrite('lfi4_predic_summarized.dat',eps1,eps2,eps3,r,ns)

! Test reheating with lnRrad and lnR

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'
  
  lnRradmin=-42
  lnRradmax = 10
  p = 2
  do i=1,npts

       lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = lfi_x_rrad(p,lnRrad,Pstar,bfoldstar)

       print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

       eps1 = lfi_epsilon_one(xstar,p)
       
!consistency test
!get lnR from lnRrad and check that it gives the same xstar
       xend = lfi_x_endinf(p)
       eps1end =  lfi_epsilon_one(xend,p)
       VendOverVstar = lfi_norm_potential(xend,p)/lfi_norm_potential(xstar,p)

       lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

       lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
       xstar = lfi_x_rreh(p,lnR,bfoldstar)
       print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

!second consistency check
!get rhoreh for chosen w and check that xstar gotten this way is the same
       w = 0._kp
       lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)
       
       xstar = lfi_x_star(p,w,lnRhoReh,Pstar,bfoldstar)
       print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
            ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar
             
    enddo

end program lfimain
