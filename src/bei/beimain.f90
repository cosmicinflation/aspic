!test the reheating derivation from slow-roll
program beimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use beisr, only : bei_epsilon_one, bei_epsilon_two,bei_epsilon_three
  use beireheat, only : bei_lnrhoend, bei_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nlambda,nbeta

  real(kp) :: lambda,beta,w,bfoldstar,lambdamin,lambdamax,betamin,betamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer


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


  do j=1,nbeta
        beta=betamin*(betamax/betamin)**(real(j,kp)/real(nbeta,kp))


    do k=1,nlambda
        lambda=lambdamin*(lambdamax/lambdamin)**(real(k,kp)/real(nlambda,kp))
     

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = bei_lnrhoend(lambda,beta,Pstar)

        print *,'lambda=',lambda,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = bei_x_star(lambda,beta,w,lnRhoReh,Pstar,bfoldstar)


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
  
        end do

     end do

 end do




end program beimain
