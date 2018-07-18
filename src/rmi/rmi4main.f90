!test the reheating derivation from slow-roll
program rmi4main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rmi4sr, only : rmi4_epsilon_one, rmi4_epsilon_two, rmi4_epsilon_three
  use rmi4reheat, only : rmi4_lnrhoreh_max, rmi4_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  use rmi4sr, only : rmi4_norm_potential, rmi4_numacc_xendmin
  use rmi4sr, only : rmi4_xendmax
  use rmi4reheat, only : rmi4_x_rreh, rmi4_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : labeps12, labnsr, labbfoldreh
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,l
  integer :: npts = 15

  integer :: Nc, Nphi0, Nxend
  real(kp) ::cmin, cmax, phi0min, phi0max, xendmin, xendmax, c, phi0, xend

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End 

  integer, parameter :: nvec = 4
  real(kp), dimension(nvec) :: cvec, phivec
  
  character(len=30), dimension(3) :: labparams
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!    Slow Roll Predictions   !!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Pstar = powerAmpScalar

  call delete_file('rmi4_predic.dat')
  call delete_file('rmi4_nsr.dat')

  labparams(1) = 'xend'
  labparams(2) = 'phi0'
  labparams(3) = 'c   '

  call aspicwrite_header('rmi4',labeps12, labnsr, labbfoldreh, labparams)

  
  
  Nxend=100

  !  w = 1._kp/3._kp
  w=0._kp


  cvec = (/-0.0005, -0.0005, -0.01, -0.01/)
  phivec = (/2.0, 10.0, 2.0, 10.0/)
  
  do k=1,nvec
     c=cvec(k)
     phi0=phivec(k)

     xendmin =  rmi4_numacc_xendmin(c,phi0) !Using an asymptotic expression for eps1 when x->1, and requiring eps1>epsilon(1._kp) for numerical convergence
     xendmax = exp(1._kp)

     do l=0,Nxend 
        xend=xendmin*(xendmax/xendmin)**(real(l,kp)/Nxend)  !logarithmic step

        if (xend.gt.rmi4_xendmax(120._kp,c,phi0)) cycle

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = rmi4_lnrhoreh_max(c,phi0,xend,Pstar)


        print *,'c=',c,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = rmi4_x_star(c,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar


           eps1 = rmi4_epsilon_one(xstar,c,phi0)
           eps2 = rmi4_epsilon_two(xstar,c,phi0)
           eps3 = rmi4_epsilon_three(xstar,c,phi0)


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)

           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('rmi4_predic.dat',c,phi0,xend,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('rmi4_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend,phi0,c/))
           
        end do

     end do

  end do

  call aspicwrite_end()
  
  ! end do

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  c = -0.001
  phi0 = 10.
  xend = 0.1_kp*rmi4_xendmax(120._kp,c,phi0)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rmi4_x_rrad(c,phi0,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rmi4_epsilon_one(xstar,c,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  rmi4_epsilon_one(xend,c,phi0)
     VendOverVstar = rmi4_norm_potential(xend,c,phi0)/rmi4_norm_potential(xstar,c,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rmi4_x_rreh(c,phi0,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rmi4_x_star(c,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program rmi4main
