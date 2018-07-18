!test the reheating derivation from slow-roll
program sdimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use sdisr, only : sdi_epsilon_one, sdi_epsilon_two, sdi_epsilon_three
  use sdireheat, only : sdi_lnrhoreh_max, sdi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use sdisr, only : sdi_norm_potential, sdi_phizeromin, sdi_efold_primitive
  use sdisr, only : sdi_numacc_xendmin, sdi_numacc_xendmax, sdi_numacc_xinimin
  use sdireheat, only : sdi_x_rreh, sdi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts, Nphi0, Nxend

  real(kp) :: phi0,xendinf,w,bfoldstar
  real(kp) :: xinimin,xendmin,xendmax,phi0min,phi0max
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), parameter :: efoldNum = 120._kp

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  integer, parameter :: nvec = 3
  real(kp), dimension(nvec) :: phivec

  Pstar = powerAmpScalar


  npts = 15

  w=0._kp
  !  w = 1._kp/3._kp

  call aspicwrite_header('sdi',labeps12,labnsr,labbfoldreh,(/'xend','phi0'/))
  
  call delete_file('sdi_predic.dat')
  call delete_file('sdi_nsr.dat')



  phivec = (/6.0,7.0,10.0/)

  Nxend = 20


  do j=1,nvec
     phi0 = phivec(j)


     xendmin = sdi_numacc_xendmin(efoldNum,phi0)
     xendmax = sdi_numacc_xendmax(phi0)
     xinimin = sdi_numacc_xinimin(phi0)

     if (phi0.lt.sdi_phizeromin()) cycle

!debug/test
     print *,'xinimin= xendmin xendmax= ',xinimin,xendmin,xendmax
     print *,'efoldMax= ',sdi_efold_primitive(xinimin,phi0) - sdi_efold_primitive(xendmin,phi0)


     xendmin=10._kp**(-3.)
     xendmax=50._kp


     do k=1,Nxend
        xendinf=xendmin*(xendmax/xendmin)**(real(k,kp)/real(Nxend,kp))

        if (xendinf.lt.sdi_numacc_xendmin(efoldNum,phi0)) cycle
        if (xendinf.gt.sdi_numacc_xendmax(phi0)) cycle

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = sdi_lnrhoreh_max(phi0,xendinf,Pstar)

        print *,'phi0=',phi0,'xendinf=',xendinf,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = sdi_x_star(phi0,xendinf,w,lnRhoReh,Pstar,bfoldstar)

           !print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

           eps1 = sdi_epsilon_one(xstar,phi0)
           eps2 = sdi_epsilon_two(xstar,phi0)
           eps3 = sdi_epsilon_three(xstar,phi0)

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           print*, 'ns=',ns,'r=',r,'bfoldstar=',bfoldstar

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xendinf,phi0/))
           
           call livewrite('sdi_predic.dat',phi0,xendinf,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('sdi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do

  call aspicwrite_end()
  
  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  phi0 = 100._kp
  xend = 20._kp
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = sdi_x_rrad(phi0,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = sdi_epsilon_one(xstar,phi0)

!consistency test
!get lnR from lnRrad and check that it gives the same xstar
     eps1end =  sdi_epsilon_one(xend,phi0)
     VendOverVstar = sdi_norm_potential(xend)/sdi_norm_potential(xstar)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = sdi_x_rreh(phi0,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

!second consistency check
!get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = sdi_x_star(phi0,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program sdimain
