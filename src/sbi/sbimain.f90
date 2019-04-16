!test the reheating derivation from slow-roll
program sbimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use sbisr, only : sbi_epsilon_one, sbi_epsilon_two, sbi_epsilon_three, sbi_alphamin
  use sbireheat, only : sbi_lnrhoreh_max, sbi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use sbisr, only : sbi_norm_potential, sbi_x_endinf
  use sbireheat, only : sbi_x_rreh, sbi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nbeta

  real(kp) :: alpha,beta,w,bfoldstar,alphamin,alphamax,betamin,betamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('sbi_predic.dat')
  call delete_file('sbi_nsr.dat')

  call aspicwrite_header('sbi',labeps12,labnsr,labbfoldreh,(/'alpha','beta '/))
  
  npts = 20

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    beta=0.00005    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  beta=0.00005
  alphamin=sbi_alphamin(beta)*(1._kp+epsilon(1._kp))
  alphamax=5000._kp*alphamin
  alphamax=60._kp*alphamin
  nalpha=100

  do k=0,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !log step

     lnRhoRehMin = lnRhoNuc
     xEnd = sbi_x_endinf(alpha,beta)
     lnRhoRehMax = sbi_lnrhoreh_max(alpha,beta,xend,Pstar)

     print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = sbi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)


        eps1 = sbi_epsilon_one(xstar,alpha,beta)
        eps2 = sbi_epsilon_two(xstar,alpha,beta)
        eps3 = sbi_epsilon_three(xstar,alpha,beta)


        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,beta/))
        
        call livewrite('sbi_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('sbi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!     beta=0.001     !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  beta=0.001
  alphamin=sbi_alphamin(beta)*(1._kp+epsilon(1._kp))
  alphamax=20._kp*alphamin
  nalpha=100


  do k=0,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !log step

     lnRhoRehMin = lnRhoNuc
     xEnd = sbi_x_endinf(alpha,beta)
     lnRhoRehMax = sbi_lnrhoreh_max(alpha,beta,xend,Pstar)

     print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=npts,1,-1

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = sbi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)


        eps1 = sbi_epsilon_one(xstar,alpha,beta)
        eps2 = sbi_epsilon_two(xstar,alpha,beta)
        eps3 = sbi_epsilon_three(xstar,alpha,beta)


        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,beta/))
        
        call livewrite('sbi_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('sbi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do

  
  call aspicwrite_end()

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha =50.*sbi_alphamin(beta)
  beta =5e-3
  xEnd = sbi_x_endinf(alpha,beta)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = sbi_x_rrad(alpha,beta,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = sbi_epsilon_one(xstar,alpha,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  sbi_epsilon_one(xend,alpha,beta)
     VendOverVstar = sbi_norm_potential(xend,alpha,beta)/sbi_norm_potential(xstar,alpha,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = sbi_x_rreh(alpha,beta,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = sbi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!     alpha=alphamin(beta)     !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('sbi_alphamin_predic.dat')
  call delete_file('sbi_alphamin_nsr.dat')

  call aspicwrite_header('sbiamin',labeps12,labnsr,labbfoldreh,(/'beta '/))
  
  betamin=0.000009
  betamax=0.1
  nbeta=100


  do k=0,nbeta
     beta=betamin*(betamax/betamin)**(real(k,kp)/real(nbeta,kp)) !log step
     alpha=sbi_alphamin(beta)

     lnRhoRehMin = lnRhoNuc
     xend = sbi_x_endinf(alpha,beta)
     lnRhoRehMax = sbi_lnrhoreh_max(alpha,beta,xend,Pstar)

     print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=npts,1,-1

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = sbi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)

        eps1 = sbi_epsilon_one(xstar,alpha,beta)
        eps2 = sbi_epsilon_two(xstar,alpha,beta)
        eps3 = sbi_epsilon_three(xstar,alpha,beta)


        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/beta/))
        
        
        call livewrite('sbi_alphamin_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('sbi_alphamin_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do

  call aspicwrite_end()
  
 write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  beta =5e-3
  alpha =50.*sbi_alphamin(beta)
  xend = sbi_x_endinf(alpha,beta)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = sbi_x_rrad(alpha,beta,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = sbi_epsilon_one(xstar,alpha,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  sbi_epsilon_one(xend,alpha,beta)
     VendOverVstar = sbi_norm_potential(xend,alpha,beta)/sbi_norm_potential(xstar,alpha,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = sbi_x_rreh(alpha,beta,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = sbi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program sbimain
