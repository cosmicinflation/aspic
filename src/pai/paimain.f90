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

  integer :: i,j,k,n
  integer :: npts, Nmu, Nxend

  real(kp) :: mu,xendinf,w,bfoldstar
  real(kp) :: xinimin,xendmin,xendmax,mumin,mumax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), parameter :: efoldNum = 120._kp

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  integer, parameter :: nvec = 3
  real(kp), dimension(nvec) :: phivec

  real(kp) :: xmin,xmax,x,V1
  
  Pstar = powerAmpScalar


  call delete_file('sdi_potential.dat')
  call delete_file('sdi_slowroll.dat')

  n=150

  xmin = 0._kp
  xmax = 7_kp
  
  do i=1,n
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)        

     V1 = sdi_norm_potential(x,mu=1._kp)


     call livewrite('sdi_potential.dat',x,V1)
     
     eps1 = sdi_epsilon_one(x,mu=1._kp)
     eps2 = sdi_epsilon_two(x,mu=1._kp)
     eps3 = sdi_epsilon_three(x,mu=1._kp)
          
     call livewrite('sdi_slowroll.dat',x,eps1,eps2,eps3)
     
  enddo



  

  npts = 15

  w=0._kp
  !  w = 1._kp/3._kp

  call aspicwrite_header('sdi',labeps12,labnsr,labbfoldreh,(/'xend','mu  '/))
  
  call delete_file('sdi_predic.dat')
  call delete_file('sdi_nsr.dat')



  phivec = (/6.0,7.0,10.0/)

  Nxend = 20


  do j=1,nvec
     mu = phivec(j)


     xendmin = sdi_numacc_xendmin(efoldNum,mu)
     xendmax = sdi_numacc_xendmax(mu)
     xinimin = sdi_numacc_xinimin(mu)

     if (mu.lt.sdi_phizeromin()) cycle

!debug/test
     print *,'xinimin= xendmin xendmax= ',xinimin,xendmin,xendmax
     print *,'efoldMax= ',sdi_efold_primitive(xinimin,mu) - sdi_efold_primitive(xendmin,mu)


     xendmin=10._kp**(-3.)
     xendmax=50._kp


     do k=1,Nxend
        xendinf=xendmin*(xendmax/xendmin)**(real(k,kp)/real(Nxend,kp))

        if (xendinf.lt.sdi_numacc_xendmin(efoldNum,mu)) cycle
        if (xendinf.gt.sdi_numacc_xendmax(mu)) cycle

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = sdi_lnrhoreh_max(mu,xendinf,Pstar)

        print *,'mu=',mu,'xendinf=',xendinf,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = sdi_x_star(mu,xendinf,w,lnRhoReh,Pstar,bfoldstar)

           !print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

           eps1 = sdi_epsilon_one(xstar,mu)
           eps2 = sdi_epsilon_two(xstar,mu)
           eps3 = sdi_epsilon_three(xstar,mu)

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           print*, 'ns=',ns,'r=',r,'bfoldstar=',bfoldstar

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xendinf,mu/))
           
           call livewrite('sdi_predic.dat',mu,xendinf,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('sdi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do

  call aspicwrite_end()
  
  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  mu = 100._kp
  xend = 20._kp
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = sdi_x_rrad(mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = sdi_epsilon_one(xstar,mu)

!consistency test
!get lnR from lnRrad and check that it gives the same xstar
     eps1end =  sdi_epsilon_one(xend,mu)
     VendOverVstar = sdi_norm_potential(xend,mu)/sdi_norm_potential(xstar,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = sdi_x_rreh(mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

!second consistency check
!get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = sdi_x_star(mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program sdimain
