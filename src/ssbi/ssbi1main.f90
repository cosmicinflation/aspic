!test the reheating derivation from slow-roll
program ssbi1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ssbi1sr, only : ssbi1_epsilon_one, ssbi1_epsilon_two, ssbi1_epsilon_three, ssbi1_alphamin
  use ssbi1reheat, only : ssbi1_lnrhoreh_max, ssbi1_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use ssbi1sr, only : ssbi1_norm_potential, ssbi1_x_endinf
  use ssbi1reheat, only : ssbi1_x_rreh, ssbi1_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 15

  integer :: Nalpha,Nbeta
  real(kp) ::alphamin, alphamax, alpha, beta, betamin, betamax

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

  call delete_file('ssbi1_alphamin.dat')

  Nbeta=1000
  betamin=0.00000001_kp!10._kp**(-3._kp)
  betamax=100._kp!10._kp**(3._kp)

  do j=0,Nbeta 
     beta=betamin*(betamax/betamin)**(real(j,kp)/Nbeta)  !logarithmic step
     !  beta=betamin+(betamax-betamin)*(real(j,kp)/Nbeta)  !arithmetic step
     alpha=ssbi1_alphamin(beta)
     call livewrite('ssbi1_alphamin.dat',beta,alpha)
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

  Nalpha=100
  !  w = 1._kp/3._kp
  w=0._kp

  Nbeta=3
  allocate(betavalues(1:Nbeta))
  betavalues(1)=10._kp**(-3._kp)
  betavalues(2)=10._kp**(-1._kp)
  betavalues(3)=10._kp**(1._kp)


  call aspicwrite_header('ssbi1',labeps12,labnsr,labbfoldreh,(/'alpha','beta '/))
  
  call delete_file('ssbi1_predic.dat')
  call delete_file('ssbi1_nsr.dat')

  do j=1,Nbeta
     beta=betavalues(j)

     alphamin=max(ssbi1_alphamin(beta)*(1._kp+100000._kp*epsilon(1._kp)),10._kp**(-3._kp))
     alphamax=alphamin*10._kp**(5.7_kp)

     do k=0,Nalpha 
        alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/Nalpha)  !logarithmic step

        lnRhoRehMin = lnRhoNuc
        xEnd = ssbi1_x_endinf(alpha,beta)
        lnRhoRehMax = ssbi1_lnrhoreh_max(alpha,beta,xend,Pstar)


        print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = ssbi1_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar


           eps1 = ssbi1_epsilon_one(xstar,alpha,beta)
           eps2 = ssbi1_epsilon_two(xstar,alpha,beta)
           eps3 = ssbi1_epsilon_three(xstar,alpha,beta)


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)

           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,beta/))
           
           call livewrite('ssbi1_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('ssbi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do

  call aspicwrite_end()

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  beta = 0.1
  alpha = 100.*ssbi1_alphamin(beta)
  xEnd = ssbi1_x_endinf(alpha,beta)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = ssbi1_x_rrad(alpha,beta,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = ssbi1_epsilon_one(xstar,alpha,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  ssbi1_epsilon_one(xend,alpha,beta)
     VendOverVstar = ssbi1_norm_potential(xend,alpha,beta)/ssbi1_norm_potential(xstar,alpha,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = ssbi1_x_rreh(alpha,beta,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = ssbi1_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program ssbi1main
