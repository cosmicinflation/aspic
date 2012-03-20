!test the reheating derivation from slow-roll
program mssmimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use mssmisr, only : mssmi_epsilon_one, mssmi_epsilon_two, mssmi_epsilon_three,mssmi_x_epsilon1_min,mssmi_alpha_min
  use mssmireheat, only : mssmi_lnrhoend, mssmi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use mssmisr, only : mssmi_efold_primitive,mssmi_x_epsilon1_min

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nbeta

  real(kp) :: alpha,beta,w,bfoldstar,alphamin,alphamax,betamin,betamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::x,xmin,xmax,Riem


  Pstar = powerAmpScalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          Calculates the prior space and               !!
!!            checks mssmi_efold_primitive               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  npts=1000
!  alphamin=0._kp
!  alphamax=1._kp
!  betamin=10._kp**(-5.)
!  betamax=1._kp

!  call delete_file('mssmi_x_eps1min.dat')
!  do i=1,npts
!       alpha=alphamin+(alphamax-alphamin)*(real(i,kp)/real(npts,kp))
!       call livewrite('mssmi_x_eps1min.dat',alpha,mssmi_x_epsilon1_min(alpha,0.1_kp), &
!                 mssmi_x_epsilon1_min(alpha,0.05_kp),mssmi_x_epsilon1_min(alpha,0.01_kp))
!  end do

!  call delete_file('mssmi_eps1min.dat')
!  do i=1,npts
!       alpha=alphamin+(alphamax-alphamin)*(real(i,kp)/real(npts,kp))
!       call livewrite('mssmi_eps1min.dat',alpha,mssmi_epsilon_one(mssmi_x_epsilon1_min(alpha,0.1_kp),alpha,0.1_kp), &
!                 mssmi_epsilon_one(mssmi_x_epsilon1_min(alpha,0.05_kp),alpha,0.05_kp), &
!                 mssmi_epsilon_one(mssmi_x_epsilon1_min(alpha,0.01_kp),alpha,0.01_kp) )
!  end do

!  call delete_file('mssmi_alphamin.dat')
!  do i=1,npts
!       beta=betamin*(betamax/betamin)**(real(i,kp)/real(npts,kp))
!       call livewrite('mssmi_alphamin.dat',beta,mssmi_alpha_min(beta))
!  end do

!  call delete_file('mssmi_PrimitiveTest.dat')
!  alpha=10._kp**(-10.)
!  beta=10._kp**(-5.)
!  xmin=0.
!  xmax=1._kp
!  Riem=0._kp
!  do i=1,npts
!       x=xmin+(xmax-xmin)*(real(i,kp)/real(npts,kp))
!       Riem=Riem+(xmax-xmin)/real(npts,kp)*((x-alpha*x**5+beta*x**9)/ &
!            (2._kp-6._kp*alpha*x**4+10._kp*beta*x**8))
!       call livewrite('mssmi_PrimitiveTest.dat',x,&
!                    mssmi_efold_primitive(x,alpha,beta)-mssmi_efold_primitive(xmin,alpha,beta), &
!                   x**2/4._kp-xmin**2/4._kp,Riem)
!  end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts = 20
  nalpha=60
  nbeta=3

  w=0._kp
!  w = 1._kp/3._kp

  call delete_file('mssmi_predic.dat')
  call delete_file('mssmi_nsr.dat')


  do j=0,nbeta
       if (j.eq.0) beta=10._kp**(-12._kp)
       if (j.eq.1) beta=3.16228*10._kp**(-11._kp)
       if (j.eq.2) beta=10._kp**(-10._kp)
       if (j.eq.3) beta=10._kp**(-9.6)

  !Prior on alpha
  alphamin=max(mssmi_alpha_min(beta),10._kp**(-10.))

  ! Value of alpha above which one cannot realize the required bfoldstar e-folds. (determined numerically)
  alphamax=sqrt(20._kp*beta/9._kp)*(1._kp-epsilon(1._kp))
  if (j.eq.0) alphamax=sqrt(20._kp*beta/9._kp)*(1._kp-epsilon(1._kp))*4._kp
  if (j.eq.1) alphamax=sqrt(20._kp*beta/9._kp)*(1._kp-epsilon(1._kp))*1.3_kp


  do k=1,nalpha
       alpha=alphamin+(alphamax-alphamin)*(real(k,kp)/real(nalpha,kp))
     

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = mssmi_lnrhoend(alpha,beta,Pstar)

  print *,'alpha=',alpha,'beta=',beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = mssmi_x_star(alpha,beta,w,lnRhoReh,Pstar,bfoldstar)


       eps1 = mssmi_epsilon_one(xstar,alpha,beta)
       eps2 = mssmi_epsilon_two(xstar,alpha,beta)
       eps3 = mssmi_epsilon_three(xstar,alpha,beta)


       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('mssmi_predic.dat',alpha,beta,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('mssmi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

  end do

 end do




end program mssmimain
