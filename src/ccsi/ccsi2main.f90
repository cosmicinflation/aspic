!test the reheating derivation from slow-roll
program ccsi2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ccsi2sr, only : ccsi2_epsilon_one, ccsi2_epsilon_two, ccsi2_epsilon_three
  use ccsi2reheat, only : ccsi2_lnrhoreh_max, ccsi2_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use ccsi2sr, only : ccsi2_norm_potential, ccsi2_numacc_xendmin
  use ccsi2reheat, only : ccsi2_x_rreh, ccsi2_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 10

  integer :: Np=3
  real(kp) :: alphamin=1d-5
  real(kp) :: alphamax=1d-3

  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad

  integer :: Ne = 10
  real(kp) :: xendmin, xendmax
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar

  call delete_file('ccsi2_predic.dat')
  call delete_file('ccsi2_nsr.dat')

  call aspicwrite_header('ccsi2',labeps12,labnsr,labbfoldreh,(/'xend ','alpha'/))

  !  w = 1._kp/3._kp
  w=0._kp


  do j=0,Np 

     alpha=10._kp**(log10(alphamin)+(log10(alphamax/alphamin))*(real(j,kp)/Np))

     print *

     print *,'alpha= ',alpha    

     lnRhoRehMin = lnRhoNuc
     xendmin = ccsi2_numacc_xendmin(120._kp,alpha)
     xendmax = 5*xendmin

     do k=0,Ne

        xend = xendmin + (xendmax-xendmin)*real(k,kp)/Ne

        lnRhoRehMax = ccsi2_lnrhoreh_max(alpha,xend,Pstar)

        print *,'xend=',xend,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = ccsi2_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

           eps1 = ccsi2_epsilon_one(xstar,alpha)
           eps2 = ccsi2_epsilon_two(xstar,alpha)
           eps3 = ccsi2_epsilon_three(xstar,alpha)


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)


           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('ccsi2_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('ccsi2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend,alpha/))
           
        end do

     enddo

  end do

  call aspicwrite_end()
  
  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 0.001
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = ccsi2_x_rrad(alpha,xend, lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = ccsi2_epsilon_one(xstar,alpha)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  ccsi2_epsilon_one(xend,alpha)
     VendOverVstar = ccsi2_norm_potential(xend,alpha)/ccsi2_norm_potential(xstar,alpha)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = ccsi2_x_rreh(alpha,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = ccsi2_x_star(alpha,xend, w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program ccsi2main
