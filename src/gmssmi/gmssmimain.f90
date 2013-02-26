!test the reheating derivation from slow-roll
program gmssmimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use gmssmisr, only : gmssmi_epsilon_one, gmssmi_epsilon_two, gmssmi_epsilon_three
  use gmssmicommon, only : gmssmi_x_epsonemin, gmssmi_x_epstwozero, gmssmi_x_epsonezero
  use gmssmisr, only : gmssmi_alphamin, gmssmi_alphamax, gmssmi_efold_primitive,gmssmi_x_endinf
  use gmssmireheat, only : gmssmi_lnrhoend, gmssmi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev


  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nphi0

  real(kp) :: alpha,phi0,w,bfoldstar,alphamin,alphamax,phi0min,phi0max
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::x,xmin,xmax,Riem
  real(kp) :: phi0A,phi0B,phi0C,phi0D,xendNumA,xendAnalA,xendNumB, &
              xendAnalB,xendNumC,xendAnalC,xendNumD,xendAnalD,alphaminAppr

  real(kp), dimension(2) :: xEpsOneZero, xEpsTwoZero

  Pstar = powerAmpScalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            checks gmssmi_efold_primitive              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  npts=1000
!  alphamin=0._kp
!  alphamax=1._kp
!  phi0min=10._kp**(-5.)
!  phi0max=1._kp

!  call delete_file('gmssmi_PrimitiveTest.dat')
!  alpha=10._kp**(-10.)
!  phi0=10._kp**(-5.)
!  xmin=0.
!  xmax=1._kp
!  Riem=0._kp
!  do i=1,npts
!       x=xmin+(xmax-xmin)*(real(i,kp)/real(npts,kp))
!       Riem=Riem+(xmax-xmin)/real(npts,kp)*((x-alpha*x**5+phi0*x**9)/ &
!            (2._kp-6._kp*alpha*x**4+10._kp*phi0*x**8))
!       call livewrite('gmssmi_PrimitiveTest.dat',x,&
!                    gmssmi_efold_primitive(x,alpha,phi0)-gmssmi_efold_primitive(xmin,alpha,phi0), &
!                   x**2/4._kp-xmin**2/4._kp,Riem)
!  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            Calculates the prior space                 !!
!!                     functions                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts=200
  alphamin=0._kp
  alphamax=2.5_kp
  phi0min=10._kp**(-6._kp)
  phi0max=10._kp**(0._kp)

  call delete_file('gmssmi_x_eps1min.dat')
  do i=1,npts
     alpha=alphamin+(alphamax-alphamin)*(real(i,kp)/real(npts,kp))
     xEpsOneZero = gmssmi_x_epsonezero(alpha)
     xEpsTwoZero = gmssmi_x_epstwozero(alpha)
     call livewrite('gmssmi_x_eps1min.dat',alpha,xEpsOneZero(2),xEpsOneZero(1) &
          , xEpsTwoZero(2), xEpsTwoZero(1),gmssmi_x_epsonemin(alpha))
  end do
  print*,'gmssmi_x_eps1min.dat written'


  call delete_file('gmssmi_eps1min.dat')
  do i=1,npts
       alpha=alphamin+(alphamax-alphamin)*(real(i,kp)/real(npts,kp))
       call livewrite('gmssmi_eps1min.dat',alpha,gmssmi_epsilon_one(gmssmi_x_epsonemin(alpha),alpha,1._kp))
  end do
  print*,'gmssmi_eps1min.dat'

 call delete_file('gmssmi_alphamin.dat')
  do i=1,npts
       phi0=phi0min*(phi0max/phi0min)**(real(i,kp)/real(npts,kp))
       alphaminAppr=1-phi0*4._kp*sqrt(2._kp)/15._kp
       call livewrite('gmssmi_alphamin.dat',phi0,gmssmi_alphamin(phi0),alphaminAppr)
  end do
  print*,'gmssmi_alphamin.dat'

  print*,'prior functions written.'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts = 20
  nphi0=4


  w=0._kp
!  w = 1._kp/3._kp

  call delete_file('gmssmi_predic.dat')
  call delete_file('gmssmi_nsr.dat')

  !Case alpha>1

  do j=0,nphi0
       
       if (j.eq.0) phi0=10._kp**(-0._kp)
       if (j.eq.1) phi0=10._kp**(-0.5_kp)
       if (j.eq.2) phi0=10._kp**(-1._kp)
       if (j.eq.3) phi0=10._kp**(-1.5_kp)
       if (j.eq.4) phi0=10._kp**(-2._kp)


       if (j.eq.0) nalpha=20
       if (j.eq.1) nalpha=20
       if (j.eq.2) nalpha=20
       if (j.eq.3) nalpha=20
       if (j.eq.4) nalpha=20

       !Prior on alpha

       alphamin=1._kp+epsilon(1._kp)
       alphamax=1._kp+2*phi0**4/60._kp**2*acos(-1._kp)**2/900._kp


       do k=0,nalpha
       alpha=alphamin+(alphamax-alphamin)*(real(k,kp)/real(nalpha,kp)) !arithmetic step
       alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !logarithmic step
     

       lnRhoRehMin = lnRhoNuc
       lnRhoRehMax = gmssmi_lnrhoend(alpha,phi0,Pstar)

       print *,'alpha=',alpha,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

       do i=1,npts

         lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

         xstar = gmssmi_x_star(alpha,phi0,w,lnRhoReh,Pstar,bfoldstar)


         eps1 = gmssmi_epsilon_one(xstar,alpha,phi0)
         eps2 = gmssmi_epsilon_two(xstar,alpha,phi0)
         eps3 = gmssmi_epsilon_three(xstar,alpha,phi0)


         print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

         logErehGeV = log_energy_reheat_ingev(lnRhoReh)
         Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

         ns = 1._kp - 2._kp*eps1 - eps2
         r =16._kp*eps1

         call livewrite('gmssmi_predic.dat',log(alpha-1._kp)/log(10._kp),sign(1._kp,alpha-1._kp),phi0,eps1,eps2,eps3,r,ns,Treh)
  
    end do

  end do

 end do

  !Case alpha<1

  do j=0,nphi0
       
       if (j.eq.0) phi0=10._kp**(-0._kp)
       if (j.eq.1) phi0=10._kp**(-0.5_kp)
       if (j.eq.2) phi0=10._kp**(-1._kp)
       if (j.eq.3) phi0=10._kp**(-1.5_kp)
       if (j.eq.4) phi0=10._kp**(-2._kp)


       if (j.eq.0) nalpha=30
       if (j.eq.1) nalpha=30
       if (j.eq.2) nalpha=30
       if (j.eq.3) nalpha=30
       if (j.eq.4) nalpha=30

       !Prior on alpha

       alphamin=1._kp-2*phi0**4/50._kp**2*acos(-1._kp)**2/900._kp
       alphamax=1._kp-epsilon(1._kp)


       do k=0,nalpha
       alpha=alphamin+(alphamax-alphamin)*(real(k,kp)/real(nalpha,kp)) !arithmetic step
       alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !logarithmic step
     

       lnRhoRehMin = lnRhoNuc
       lnRhoRehMax = gmssmi_lnrhoend(alpha,phi0,Pstar)

       print *,'alpha=',alpha,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

       do i=1,npts

         lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

         xstar = gmssmi_x_star(alpha,phi0,w,lnRhoReh,Pstar,bfoldstar)


         eps1 = gmssmi_epsilon_one(xstar,alpha,phi0)
         eps2 = gmssmi_epsilon_two(xstar,alpha,phi0)
         eps3 = gmssmi_epsilon_three(xstar,alpha,phi0)


         print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

         logErehGeV = log_energy_reheat_ingev(lnRhoReh)
         Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

         ns = 1._kp - 2._kp*eps1 - eps2
         r =16._kp*eps1

         call livewrite('gmssmi_predic.dat',log(1._kp-alpha)/log(10._kp),sign(1._kp,alpha-1._kp),phi0,eps1,eps2,eps3,r,ns,Treh)
  
    end do

  end do

 end do




end program gmssmimain
