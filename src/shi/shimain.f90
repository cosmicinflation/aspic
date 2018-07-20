!test the reheating derivation from slow-roll
program shimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use shisr, only : shi_epsilon_one, shi_epsilon_two, shi_epsilon_three
  use shireheat, only : shi_lnrhoreh_max, shi_x_star
  use infinout, only : delete_file, livewrite

  use srreheat, only : log_energy_reheat_ingev 
  use specialinf, only : lambert 

  use shisr, only : shi_norm_potential, shi_x_endinf, shi_efold_primitive
  use shireheat, only : shi_x_rreh, shi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : labeps12, labnsr, labbfoldreh
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 20

  integer :: Nphi0
  real(kp) :: phi0min,phi0max

  integer :: Nalpha
  real(kp) :: alphamin,alphamax

  real(kp) :: alpha,phi0,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB

  real(kp) ::xend,xendapprox

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End

  integer, parameter :: nvec = 4
  real(kp), dimension(nvec) :: phivec
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      Reheating Predictions       !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
  Pstar = powerAmpScalar

  call delete_file('shi_predic.dat')

  call aspicwrite_header('shi',labeps12,labnsr,labbfoldreh,(/'alpha','phi0 '/))

  !  w = 1._kp/3._kp
  w=0._kp



  alphamin = 0.1
  alphamax = 20.
  Nalpha = 500
  
  phivec = (/10.0, 15.0, 20.0, 25.0/)
  
  do j=1,nvec
     phi0 = phivec(j)

     do k=0,Nalpha
        alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/Nalpha) !logarithmic step

        lnRhoRehMin = lnRhoNuc
        xEnd = shi_x_endinf(alpha,phi0)
        lnRhoRehMax = shi_lnrhoreh_max(alpha,phi0,xend,Pstar)

        print *,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = shi_x_star(alpha,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

           eps1 = shi_epsilon_one(xstar,alpha,phi0)
           eps2 = shi_epsilon_two(xstar,alpha,phi0)
           eps3 = shi_epsilon_three(xstar,alpha,phi0)

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('shi_predic.dat',alpha,phi0,eps1,eps2,eps3,r,ns,Treh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,phi0/))

        end do

     end do

  enddo

  call aspicwrite_end()
  


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 1.
  phi0 = 10.
  xEnd = shi_x_endinf(alpha,phi0)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = shi_x_rrad(alpha,phi0,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = shi_epsilon_one(xstar,alpha,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  shi_epsilon_one(xend,alpha,phi0)
     VendOverVstar = shi_norm_potential(xend,alpha,phi0)/shi_norm_potential(xstar,alpha,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = shi_x_rreh(alpha,phi0,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = shi_x_star(alpha,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program shimain
