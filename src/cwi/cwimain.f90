!test the reheating derivation from slow-roll
program cwimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use cwisr, only : cwi_epsilon_one, cwi_epsilon_two, cwi_epsilon_three, cwi_x_endinf
  use cwireheat, only : cwi_lnrhoend, cwi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev 
  use specialinf, only : lambert 

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  integer :: NQ
  real(kp) :: Qmin,Qmax

  real(kp) :: alpha,Q,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB

  real(kp) ::xend,xendapprox

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!      Tests the approximated      !!!!
  !!!!         formula for xend         !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   Qmin=0.000005
   Qmax=0.005

   call delete_file('cwi_xend.dat')

   NQ=100
   alpha=4._kp*(exp(1._kp))

   do j=0,NQ 
   Q=Qmin*(Qmax/Qmin)**(real(j,kp)/NQ) !logarithmic step
   xend=cwi_x_endinf(alpha,Q)
   xendapprox=exp(-0.25_kp)*exp(lambert(-3._kp*sqrt(2._kp)/(4._kp*alpha)*Q*exp(0.75_kp),-1)/3._kp)
   call livewrite('cwi_xend.dat',Q,xend,xendapprox)
   end do


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!      Reheating Predictions       !!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  Qmin=1
  Qmax=100.

  !!!!!!!!!!!!!!!!!!!!!
  !! Physical Regime !!
  !!!!!!!!!!!!!!!!!!!!!

  Qmin=0.00001
  Qmax=0.001

  NQ=15
  alpha=4._kp*exp(1._kp)

  Pstar = powerAmpScalar

  call delete_file('cwi_predic1.dat')

!  w = 1._kp/3._kp
  w=0._kp

  do j=0,NQ 
  Q=Qmin*(Qmax/Qmin)**(real(j,kp)/NQ) !logarithmic step

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = cwi_lnrhoend(alpha,Q,Pstar)

  print *,'Q=',Q,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)
       
       xstar = cwi_x_star(alpha,Q,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar 

       eps1 = cwi_epsilon_one(xstar,alpha,Q)
       eps2 = cwi_epsilon_two(xstar,alpha,Q)
       eps3 = cwi_epsilon_three(xstar,alpha,Q)

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('cwi_predic1.dat',alpha,Q,eps1,eps2,eps3,r,ns,Treh)
  
    end do

 end do

  !!!!!!!!!!!!!!!!!!
  !! Other Regime !!
  !!!!!!!!!!!!!!!!!!

  Qmin=1.
  Qmax=400.

  NQ=30
  alpha=4._kp*exp(1._kp)

  Pstar = powerAmpScalar

  call delete_file('cwi_predic2.dat')

!  w = 1._kp/3._kp
  w=0._kp

  do j=0,NQ 
  Q=Qmin*(Qmax/Qmin)**(real(j,kp)/NQ) !logarithmic step

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = cwi_lnrhoend(alpha,Q,Pstar)

  print *,'Q=',Q,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)
       
       xstar = cwi_x_star(alpha,Q,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar 

       eps1 = cwi_epsilon_one(xstar,alpha,Q)
       eps2 = cwi_epsilon_two(xstar,alpha,Q)
       eps3 = cwi_epsilon_three(xstar,alpha,Q)

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('cwi_predic2.dat',alpha,Q,eps1,eps2,eps3,r,ns,Treh)
  
    end do

 end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('cwi_predic_summarized.dat') 
         NQ=1000
         Qmin=1.
         Qmax=100.
         w=0._kp
         Pstar = powerAmpScalar
         do j=1,NQ
         Q=Qmin*(Qmax/Qmin)**(real(j,kp)/real(NQ,kp))
         alpha=4._kp*(exp(1._kp))
         lnRhoReh = lnRhoNuc
         xstarA = cwi_x_star(alpha,Q,w,lnRhoReh,Pstar,bfoldstar)
         eps1A = cwi_epsilon_one(xstarA,alpha,Q)
         eps2A = cwi_epsilon_two(xstarA,alpha,Q)
         eps3A = cwi_epsilon_three(xstarA,alpha,Q)
         nsA = 1._kp - 2._kp*eps1A - eps2A
         rA = 16._kp*eps1A
         lnRhoReh = cwi_lnrhoend(alpha,Q,Pstar)
         xstarB = cwi_x_star(alpha,Q,w,lnRhoReh,Pstar,bfoldstar)
         eps1B = cwi_epsilon_one(xstarB,alpha,Q)
         eps2B = cwi_epsilon_two(xstarB,alpha,Q)
         eps3B = cwi_epsilon_three(xstarB,alpha,Q)
         nsB = 1._kp - 2._kp*eps1B - eps2B
         rB =16._kp*eps1B
         call livewrite('cwi_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
         enddo


end program cwimain
