!test the reheating derivation from slow-roll
program deimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use deisr, only : dei_epsilon_one, dei_epsilon_two, dei_epsilon_three
  use deireheat, only : dei_lnrhoreh_max, dei_x_star
  use infinout, only : delete_file, livewrite

  use srreheat, only : log_energy_reheat_ingev 
  use specialinf, only : lambert 

  use deisr, only : dei_norm_potential, dei_x_endinf, dei_efold_primitive
  use deireheat, only : dei_x_rreh, dei_x_rrad
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

  integer :: Nbeta
  real(kp) :: betamin,betamax

  real(kp) :: beta,phi0,w,bfoldstar
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

  call delete_file('dei_predic.dat')

  call aspicwrite_header('dei',labeps12,labnsr,labbfoldreh,(/'beta','phi0'/))

  !  w = 1._kp/3._kp
  w=0._kp



  betamin = 0.01
  betamax = 0.99
  Nbeta = 100
  
  phivec = (/10.0, 20.0, 50.0, 100.0/)
  
  do j=1,nvec
     phi0 = phivec(j)

     do k=0,Nbeta
        beta=betamin+(betamax-betamin)*(real(k,kp)/Nbeta) !arithmetic step

        lnRhoRehMin = lnRhoNuc
        xEnd = dei_x_endinf(beta,phi0)
        lnRhoRehMax = dei_lnrhoreh_max(beta,phi0,xend,Pstar)

        print *,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = dei_x_star(beta,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

           eps1 = dei_epsilon_one(xstar,beta,phi0)
           eps2 = dei_epsilon_two(xstar,beta,phi0)
           eps3 = dei_epsilon_three(xstar,beta,phi0)

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('dei_predic.dat',beta,phi0,eps1,eps2,eps3,r,ns,Treh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/beta,phi0/))

        end do

     end do

  enddo

  call aspicwrite_end()
  


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  beta = 0.5
  phi0 = 50.
  xEnd = dei_x_endinf(beta,phi0)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = dei_x_rrad(beta,phi0,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = dei_epsilon_one(xstar,beta,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  dei_epsilon_one(xend,beta,phi0)
     VendOverVstar = dei_norm_potential(xend,beta,phi0)/dei_norm_potential(xstar,beta,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = dei_x_rreh(beta,phi0,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = dei_x_star(beta,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program deimain
