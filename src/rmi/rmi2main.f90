!test the reheating derivation from slow-roll
program rmi2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rmi2sr, only : rmi2_epsilon_one, rmi2_epsilon_two, rmi2_epsilon_three, rmi2_numacc_xendmin
  use rmi2reheat, only : rmi2_lnrhoreh_max, rmi2_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use rmi2sr, only : rmi2_norm_potential
  use rmi2reheat, only : rmi2_x_rreh, rmi2_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : labeps12, labnsr, labbfoldreh
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,l
  integer :: npts = 20

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

  call delete_file('rmi2_predic.dat')
  call delete_file('rmi2_nsr.dat')

  labparams(1) = 'xend'
  labparams(2) = 'phi0'
  labparams(3) = 'c   '
  
  call aspicwrite_header('rmi2',labeps12, labnsr, labbfoldreh, labparams)


  Nxend=100

  !  w = 1._kp/3._kp
  w=0._kp


  cvec = (/0.005, 0.005, 0.01, 0.01/)
  phivec = (/3.0, 10.0, 5.0, 10.0/)

  
  do k=1,nvec
     c=cvec(k)
     phi0=phivec(k)
  
     xendmin = rmi2_numacc_xendmin(70._kp,c,phi0)
     xendmax = exp(1._kp)


     if (xendmax .lt. xendmin) then
        print*,'xendmax<xendmin !!: not a sufficient number of efold can be realized in the region where the potential is valid!'
        cycle
     endif



     do l=0,Nxend 
        !xend=xendmin+(xendmax-xendmin)*(real(l,kp)/Nxend)  !arithmetic step
        xend=xendmin*(xendmax/xendmin)**(real(l,kp)/Nxend)  !logarithmic step
 !ultralogarithmic step
        !xend=exp(exp((real(l,kp)/Nxend)*log(log(xendmax)/log(xendmin)))*log(xendmin))
 !tangent step
        ! xend=xendmin+(xendmax-xendmin)*atan(real(l,kp)/Nxend*5._kp)*2._kp/acos(-1._kp)
        
        if (rmi2_epsilon_one(xend,c,phi0).gt.1._kp) then
           print *,'c too large, eps1>1',xend,c,phi0
           cycle
        endif

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = rmi2_lnrhoreh_max(c,phi0,xend,Pstar)


        print *,'c=',c,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax, 'xend=',xend

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = rmi2_x_star(c,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = rmi2_epsilon_one(xstar,c,phi0)
           eps2 = rmi2_epsilon_two(xstar,c,phi0)
           eps3 = rmi2_epsilon_three(xstar,c,phi0)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1,'eps2star=',eps2


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)

           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('rmi2_predic.dat',c,phi0,xend,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('rmi2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

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
  c = 0.001
  phi0 = 10.
  xend = rmi2_numacc_xendmin(70._kp,c,phi0)+1.
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rmi2_x_rrad(c,phi0,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rmi2_epsilon_one(xstar,c,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  rmi2_epsilon_one(xend,c,phi0)
     VendOverVstar = rmi2_norm_potential(xend,c,phi0)/rmi2_norm_potential(xstar,c,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rmi2_x_rreh(c,phi0,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rmi2_x_star(c,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program rmi2main
