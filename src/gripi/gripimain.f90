!test the reheating derivation from slow-roll
program gripimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use gripisr, only : gripi_epsilon_one, gripi_epsilon_two, gripi_epsilon_three
  use gripisr, only : gripi_x_epsonemin, gripi_x_epstwozero, gripi_x_epsonezero
  use gripisr, only : gripi_alphamin, gripi_alphamax, gripi_efold_primitive,gripi_x_endinf
  use gripireheat, only : gripi_lnrhoreh_max, gripi_x_star
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

  real(kp), dimension(2) :: xEpsOneZero,xEpsTwoZero

  real(kp) :: alphaminAppr

  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!            Calculates the prior space                 !!
!!                     functions                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts=1000
  alphamin=0._kp
  alphamax=2.5_kp
  phi0min=10._kp**(-6._kp)
  phi0max=10._kp**(0._kp)

  call delete_file('gripi_x_eps1min.dat')
  do i=1,npts
     alpha=alphamin+(alphamax-alphamin)*(real(i,kp)/real(npts,kp))
     xEpsOneZero = gripi_x_epsonezero(alpha)
     xEpsTwoZero = gripi_x_epstwozero(alpha)

     call livewrite('gripi_x_eps1min.dat',alpha,xEpsOneZero(2),xEpsOneZero(1), &
        xEpsTwoZero(1),xEpsTwoZero(2),gripi_x_epsonemin(alpha))
  end do
  print*,'gripi_x_eps1min.dat written'


  call delete_file('gripi_eps1min.dat')
  do i=1,npts
       alpha=alphamin+(alphamax-alphamin)*(real(i,kp)/real(npts,kp))
       call livewrite('gripi_eps1min.dat',alpha,gripi_epsilon_one(gripi_x_epsonemin(alpha),alpha,1._kp))
  end do
  print*,'gripi_eps1min.dat'

! call delete_file('gripi_alphamin.dat')
!  do i=1,npts
!       phi0=phi0min*(phi0max/phi0min)**(real(i,kp)/real(npts,kp))
!       alphaminAppr=1-phi0*4._kp*sqrt(2._kp)/15._kp
!       call livewrite('gripi_alphamin.dat',phi0,gripi_alphamin(phi0),alphaminAppr)
!  end do
!  print*,'gripi_alphamin.dat'

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

  call delete_file('gripi_predic.dat')
  call delete_file('gripi_nsr.dat')

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
       alphamax=1._kp+2*phi0**4/60._kp**2*acos(-1._kp)**2/576._kp


       do k=0,nalpha
       alpha=alphamin+(alphamax-alphamin)*(real(k,kp)/real(nalpha,kp)) !arithmetic step
       alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !logarithmic step
!       alpha=1._kp-(1._kp-alphamin)*((1._kp-alphamax)/(1._kp-alphamin))** &
!            (real(k,kp)/real(nalpha,kp))!logarithmic step on 1-alpha
     

       lnRhoRehMin = lnRhoNuc
       lnRhoRehMax = gripi_lnrhoreh_max(alpha,phi0,Pstar)

       print *,'alpha=',alpha,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

       do i=1,npts

         lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

         xstar = gripi_x_star(alpha,phi0,w,lnRhoReh,Pstar,bfoldstar)


         eps1 = gripi_epsilon_one(xstar,alpha,phi0)
         eps2 = gripi_epsilon_two(xstar,alpha,phi0)
         eps3 = gripi_epsilon_three(xstar,alpha,phi0)


         print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

         logErehGeV = log_energy_reheat_ingev(lnRhoReh)
         Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

         ns = 1._kp - 2._kp*eps1 - eps2
         r =16._kp*eps1

         call livewrite('gripi_predic.dat',log(alpha-1._kp)/log(10._kp),sign(1._kp,alpha-1._kp),phi0,eps1,eps2,eps3,r,ns,Treh)
  
    end do

  end do

 end do

print*,'case alpha>1 done'


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

       alphamin=1._kp-2*phi0**4/50._kp**2*acos(-1._kp)**2/576._kp
       alphamax=1._kp-10.*epsilon(1._kp)


       do k=0,nalpha
       alpha=alphamin+(alphamax-alphamin)*(real(k,kp)/real(nalpha,kp)) !arithmetic step
       alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !logarithmic step
!       alpha=1._kp-(1._kp-alphamin)*((1._kp-alphamax)/(1._kp-alphamin))** &
!            (real(k,kp)/real(nalpha,kp))!logarithmic step on 1-alpha
     

       lnRhoRehMin = lnRhoNuc
       lnRhoRehMax = gripi_lnrhoreh_max(alpha,phi0,Pstar)

       print *,'alpha=',alpha,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

       do i=1,npts

         lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

         xstar = gripi_x_star(alpha,phi0,w,lnRhoReh,Pstar,bfoldstar)


         eps1 = gripi_epsilon_one(xstar,alpha,phi0)
         eps2 = gripi_epsilon_two(xstar,alpha,phi0)
         eps3 = gripi_epsilon_three(xstar,alpha,phi0)


         print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

         logErehGeV = log_energy_reheat_ingev(lnRhoReh)
         Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

         ns = 1._kp - 2._kp*eps1 - eps2
         r =16._kp*eps1

         call livewrite('gripi_predic.dat',log(1._kp-alpha)/log(10._kp),sign(1._kp,alpha-1._kp),phi0,eps1,eps2,eps3,r,ns,Treh)
  
    end do

  end do

 end do

print*,'case alpha<1 done'



end program gripimain