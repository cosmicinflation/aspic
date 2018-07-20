!test the reheating derivation from slow-roll
program ssbi6main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ssbi6sr, only : ssbi6_epsilon_one, ssbi6_epsilon_two, ssbi6_epsilon_three, ssbi6_alphamax
  use ssbi6reheat, only : ssbi6_lnrhoreh_max, ssbi6_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use ssbi6sr, only : ssbi6_norm_potential, ssbi6_x_endinf
  use ssbi6reheat, only : ssbi6_x_rreh, ssbi6_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 15

  integer :: Nalpha,Nbeta
  real(kp) ::alphamin, alphamax, betamin, betamax, alpha, beta

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer
  real(kp), dimension(:), allocatable :: betavalues

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!         Prior Space        !!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('ssbi6_abs_alpha_min.dat')

  Nbeta=10000
  betamin=0.000001_kp
  betamax=100._kp

  do j=0,Nbeta 
     beta=betamin*(betamax/betamin)**(real(j,kp)/Nbeta)  !logarithmic step
     !  beta=betamin+(betamax-betamin)*(real(j,kp)/Nbeta)  !arithmetic step
     alpha=ssbi6_alphamax(beta)
     call livewrite('ssbi6_alphamax.dat',beta,alpha)
  end do

  print *,'Priors Written'





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!    Slow Roll Predictions   !!!!!!!!!!
!!!!!!!                            !!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Pstar = powerAmpScalar

  call aspicwrite_header('ssbi6',labeps12,labnsr,labbfoldreh,(/'alpha','beta '/))
  
  call delete_file('ssbi6_predic.dat')
  call delete_file('ssbi6_nsr.dat')

  Nalpha=100

  !  w = 1._kp/3._kp
  w=0._kp

  Nbeta=3
  allocate(betavalues(1:Nbeta))
  betavalues(1)=10._kp**(-5._kp)
  betavalues(2)=0.1_kp
  betavalues(3)=10._kp**(0._kp)

  do j=1,Nbeta
     beta=betavalues(j)

     alphamin=ssbi6_alphamax(beta)*1.001_kp
     if (alphamin .eq. 0._kp) alphamin=-10._kp**(-4._kp)
     alphamax=alphamin*10._kp**(4._kp)
     if (beta .eq. 10._kp**(-6._kp))   alphamax=alphamin*10._kp**(1.5_kp)

     do k=0,Nalpha 
        alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/Nalpha)  !logarithmic step

        lnRhoRehMin = lnRhoNuc
        xEnd = ssbi6_x_endinf(alpha,beta)
        lnRhoRehMax = ssbi6_lnrhoreh_max(alpha,beta,xend,Pstar)


        print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = ssbi6_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar


           eps1 = ssbi6_epsilon_one(xstar,alpha,beta)
           eps2 = ssbi6_epsilon_two(xstar,alpha,beta)
           eps3 = ssbi6_epsilon_three(xstar,alpha,beta)


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)

           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,beta/))
           
           call livewrite('ssbi6_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('ssbi6_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do

  call aspicwrite_end()
  
  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  beta = 0.1
  alpha = -0.9
  xEnd = ssbi6_x_endinf(alpha,beta)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = ssbi6_x_rrad(alpha,beta,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = ssbi6_epsilon_one(xstar,alpha,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  ssbi6_epsilon_one(xend,alpha,beta)
     VendOverVstar = ssbi6_norm_potential(xend,alpha,beta)/ssbi6_norm_potential(xstar,alpha,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = ssbi6_x_rreh(alpha,beta,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = ssbi6_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program ssbi6main
