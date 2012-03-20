!test the reheating derivation from slow-roll
program psnimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use psnisr, only : psni_epsilon_one, psni_epsilon_two, psni_epsilon_three
  use psnireheat, only : psni_lnrhoend, psni_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev


  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,k
  integer :: npts,nalpha

  real(kp) :: alpha,mu,w,bfoldstar,alphamin,alphamax,alphaminlog,alphamaxlog
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
  nalpha=20

  mu=2._kp
  mu=1._kp
  mu=0.5_kp

  w=0._kp
!  w = 1._kp/3._kp

  call delete_file('psni_predic.dat')
  call delete_file('psni_nsr.dat')


      
       alphamax=10._kp**(-1._kp)
       alphamin=10._kp**(-3._kp)

       alphaminlog=10._kp**(-7._kp)
       alphamaxlog=10._kp**(-2.3_kp)

  do k=1,2*nalpha
       if(k.lt.nalpha) then
       alpha=alphaminlog*(alphamaxlog/alphaminlog)**(real(k,kp)/real(nalpha,kp)) !log step between alphaminlog and alphamaxlog
       else
       alpha=alphamin+(alphamax-alphamin)*(real(k-nalpha,kp)/real(nalpha,kp)) !arithmetic step between alphamin and alphamax
       endif
     

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = psni_lnrhoend(alpha,mu,Pstar)

  print *,'alpha=',alpha,'mu=',mu,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = psni_x_star(alpha,mu,w,lnRhoReh,Pstar,bfoldstar)


       eps1 = psni_epsilon_one(xstar,alpha,mu)
       eps2 = psni_epsilon_two(xstar,alpha,mu)
       eps3 = psni_epsilon_three(xstar,alpha,mu)


       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('psni_predic.dat',alpha,mu,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('psni_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do



 end do




end program psnimain
