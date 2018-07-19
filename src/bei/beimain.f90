!test the reheating derivation from slow-roll
program beimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use beisr, only : bei_epsilon_one, bei_epsilon_two,bei_epsilon_three
  use beireheat, only : bei_lnrhoreh_max, bei_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use beisr, only : bei_norm_potential, bei_x_endinf
  use beireheat, only : bei_x_rreh, bei_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nlambda,nbeta

  real(kp) :: lambda,beta,w,bfoldstar,lambdamin,lambdamax,betamin,betamax
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

  npts = 20
  nlambda=3
  nbeta=10



  lambdamin=10._kp**(-6._kp)
  lambdamax=10._kp**(3._kp)

  betamin=10._kp**(-2._kp)
  betamax=10._kp**(2.5_kp)

  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('bei_predic.dat')
  call delete_file('bei_nsr.dat')

  call aspicwrite_header('bei',labeps12,labnsr,labbfoldreh,(/'beta  ','lambda'/))

  do k=1,nlambda
     lambda=lambdamin*(lambdamax/lambdamin)**(real(k,kp)/real(nlambda,kp))
  
  
     do j=1,nbeta
        beta=betamin*(betamax/betamin)**(real(j,kp)/real(nbeta,kp))
             

        lnRhoRehMin = lnRhoNuc
        xEnd=bei_x_endinf(lambda,beta) 
        lnRhoRehMax = bei_lnrhoreh_max(lambda,beta,xend,Pstar)

        print *,'lambda=',lambda,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = bei_x_star(lambda,beta,xend,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = bei_epsilon_one(xstar,lambda,beta)
           eps2 = bei_epsilon_two(xstar,lambda,beta)
           eps3 = bei_epsilon_three(xstar,lambda,beta)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1, &
                '  xstaranal=',1._kp/(lambda*beta)-sqrt(1/(2._kp*beta**2)-2._kp*bfoldstar/beta), &
                '   epsanal=',1._kp/(1._kp-4._kp*beta*bfoldstar)


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('bei_predic.dat',lambda,beta,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('bei_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/beta,lambda/))
           
        end do

     end do

  end do

  call aspicwrite_end()
  

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  lambda = 1
  beta = 2
  xEnd=bei_x_endinf(lambda,beta)
  
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = bei_x_rrad(lambda,beta,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = bei_epsilon_one(xstar,lambda,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar

     eps1end =  bei_epsilon_one(xend,lambda,beta)
     VendOverVstar = bei_norm_potential(xend,lambda,beta)/bei_norm_potential(xstar,lambda,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = bei_x_rreh(lambda,beta,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = bei_x_star(lambda,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo




end program beimain
