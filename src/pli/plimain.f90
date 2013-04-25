!test the reheating derivation from slow-roll
program plimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use plisr, only : pli_epsilon_one, pli_epsilon_two, pli_epsilon_three
  use plireheat, only : pli_lnrhoreh_max, pli_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i
  integer :: npts = 2

  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp)  ::alphamin,alphamax,eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  integer :: nalpha

  Pstar = powerAmpScalar

  call delete_file('pli_predic.dat')
  call delete_file('pli_nsr.dat')

  alpha=sqrt(2._kp)*10**(-5./2.)
  do while (alpha.le.0.9_kp)

  alpha = alpha*1.1
!  w = 1._kp/3._kp
  w=0._kp
 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = pli_lnrhoreh_max(alpha,Pstar)

  print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = pli_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

       eps1 = pli_epsilon_one(xstar,alpha)
       eps2 = pli_epsilon_two(xstar,alpha)
       eps3 = pli_epsilon_three(xstar,alpha)

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('pli_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('pli_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

  enddo

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('pli_predic_summarized.dat') 
         nalpha=1000
         alphamin=10._kp**(-3.)
         alphamax=10._kp**(0.)
         w=0._kp
         do i=1,nalpha
         alpha=alphamin*(alphamax/alphamin)**(real(i,kp)/real(nalpha,kp))
         lnRhoReh = lnRhoNuc
         xstarA = pli_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
         eps1A = pli_epsilon_one(xstarA,alpha)
         eps2A = pli_epsilon_two(xstarA,alpha)
         eps3A = pli_epsilon_three(xstarA,alpha)
         nsA = 1._kp - 2._kp*eps1A - eps2A
         rA = 16._kp*eps1A
         lnRhoReh = pli_lnrhoreh_max(alpha,Pstar)
         xstarB = pli_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
         eps1B = pli_epsilon_one(xstarB,alpha)
         eps2B = pli_epsilon_two(xstarB,alpha)
         eps3B = pli_epsilon_three(xstarB,alpha)
         nsB = 1._kp - 2._kp*eps1B - eps2B
         rB =16._kp*eps1B
         call livewrite('pli_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
         enddo



end program plimain
