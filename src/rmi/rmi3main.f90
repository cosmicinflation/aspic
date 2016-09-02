!test the reheating derivation from slow-roll
program rmi3main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rmi3sr, only : rmi3_epsilon_one, rmi3_epsilon_two, rmi3_epsilon_three
  use rmi3reheat, only : rmi3_lnrhoreh_max, rmi3_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  use infinout, only : labeps12, labnsr, labbfoldreh
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use rmi3sr, only : rmi3_norm_potential, rmi3_numacc_xendmax, rmi3_numacc_xendmin
  use rmi3reheat, only : rmi3_x_rreh, rmi3_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,l
  integer :: npts = 4

  integer :: Nc, Nphi0, Nxend
  real(kp) ::cmin, cmax, phi0min, phi0max, xendmin, xendmax, c, phi0, xend

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End 

  character(len=30), dimension(3) :: labparams
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!    Slow Roll Predictions   !!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Pstar = powerAmpScalar

  call delete_file('rmi3_predic.dat')
  call delete_file('rmi3_nsr.dat')
  
  labparams(1) = 'xend'
  labparams(2) = 'phi0'
  labparams(3) = 'c   '
  
  call aspicwrite_header('rmi3',labeps12, labnsr, labbfoldreh, labparams)
  
  Nc=10
  Nphi0=8
  Nxend=50

  !  w = 1._kp/3._kp
  w=0._kp

!  cmin=10._kp**(-3._kp)
!  cmax=10._kp**(0._kp)

  ! do j=0,Nc 
  ! c=cmin*(cmax/cmin)**(real(j,kp)/Nc)  !logarithmic step

  c=-10._kp**(-2._kp)
  !  c=-5.*10._kp**(-3._kp)
  !  c=-10._kp**(-2._kp)

  phi0max=1._kp/sqrt(abs(c))
  phi0min=phi0max/(50.)

  do k=0,Nphi0 
     phi0=phi0min*(phi0max/phi0min)**(real(k,kp)/Nphi0)  !logarithmic step

     xendmin = exp(-1._kp)
     xendmax = rmi3_numacc_xendmax(c,phi0) !Using an asymptotic expression for eps1 when x->1, and requiring eps1>epsilon(1._kp) for numerical convergence

     do l=0,Nxend 
        xend=xendmin*(xendmax/xendmin)**(real(l,kp)/Nxend)  !logarithmic step

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = rmi3_lnrhoreh_max(c,phi0,xend,Pstar)


        print *,'c=',c,'phi0=',phi0,'xend=', xend, 'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = rmi3_x_star(c,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = rmi3_epsilon_one(xstar,c,phi0)
           eps2 = rmi3_epsilon_two(xstar,c,phi0)
           eps3 = rmi3_epsilon_three(xstar,c,phi0)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1,'eps2star=',eps2


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)

           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('rmi3_predic.dat',c,phi0,xend,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('rmi3_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

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
  xend = 10*exp(-1._kp)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rmi3_x_rrad(c,phi0,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rmi3_epsilon_one(xstar,c,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  rmi3_epsilon_one(xend,c,phi0)
     VendOverVstar = rmi3_norm_potential(xend,c,phi0)/rmi3_norm_potential(xstar,c,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rmi3_x_rreh(c,phi0,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rmi3_x_star(c,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program rmi3main
