!test the reheating derivation from slow-roll
program mssmimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use mssmisr, only : mssmi_epsilon_one, mssmi_epsilon_two, mssmi_epsilon_three, mssmi_x_epsonemin
  use mssmireheat, only : mssmi_lnrhoend, mssmi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev


  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts,nalpha

  real(kp) :: alpha,w,bfoldstar,alphamin,alphamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::x,xmin,xmax,Riem


  real(kp) :: eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB


  Pstar = powerAmpScalar



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts = 20
  nalpha=30

  w=0._kp
!  w = 1._kp/3._kp

  call delete_file('mssmi_predic.dat')
  call delete_file('mssmi_nsr.dat')


  !Prior on alpha
  alphamin=10._kp**(-8.)
  alphamax=10._kp**(-3.)


  do j=0,nalpha
       alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))
 

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = mssmi_lnrhoend(alpha,Pstar)

  print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = mssmi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)


       eps1 = mssmi_epsilon_one(xstar,alpha)
       eps2 = mssmi_epsilon_two(xstar,alpha)
       eps3 = mssmi_epsilon_three(xstar,alpha)


       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('mssmi_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('mssmi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do


 end do

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('mssmi_predic_summarized.dat') 
         nalpha=1000
         alphamin=10._kp**(-8.)
         alphamax=10._kp**(-3.)
         w=0._kp
         do j=1,nalpha
         alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))
         lnRhoReh = lnRhoNuc
         xstarA = mssmi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
         eps1A = mssmi_epsilon_one(xstarA,alpha)
         eps2A = mssmi_epsilon_two(xstarA,alpha)
         eps3A = mssmi_epsilon_three(xstarA,alpha)
         nsA = 1._kp - 2._kp*eps1A - eps2A
         rA = 16._kp*eps1A
         lnRhoReh = mssmi_lnrhoend(alpha,Pstar)
         xstarB = mssmi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
         eps1B = mssmi_epsilon_one(xstarB,alpha)
         eps2B = mssmi_epsilon_two(xstarB,alpha)
         eps3B = mssmi_epsilon_three(xstarB,alpha)
         nsB = 1._kp - 2._kp*eps1B - eps2B
         rB =16._kp*eps1B
         if ((rA .gt. 0._kp) .and. (nsA .gt.0._kp) .and. &
                  (rB .gt. 0._kp) .and. (nsB .gt. 0._kp) &
                  .and. (eps1A .gt. 10.**(-10.)) .and. (eps1B .gt. 10.**(-10.)) &
                 .and. (eps2A .gt. 10.**(-10.)) .and. (eps2B .gt. 10.**(-10.))) then !to remove NaN and error points
         call livewrite('mssmi_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
         end if
         enddo


end program mssmimain
