!test the reheating derivation from slow-roll
program limain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lisr, only : li_epsilon_one, li_epsilon_two, li_epsilon_three
  use lireheat, only : li_lnrhoreh_max, li_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts

  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer
  
  real(kp) ::alphamin,alphamax

  real(kp) :: eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  integer :: nalpha

  Pstar = powerAmpScalar

  call delete_file('li_predic.dat')
  call delete_file('li_nsr.dat')

!!!!!!!!!!!!!!!!!!
!!!  alpha>0   !!!
!!!!!!!!!!!!!!!!!!

 npts = 20

  alphamin=0.002
  alphamax=100000._kp
  nalpha=20

!  w = 1._kp/3._kp
  w=0._kp

    do k=0,nalpha
      alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp))

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = li_lnrhoreh_max(alpha,Pstar)

  print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = li_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

       eps1 = li_epsilon_one(xstar,alpha)
       eps2 = li_epsilon_two(xstar,alpha)
       eps3 = li_epsilon_three(xstar,alpha)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1*=',eps1

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('li_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('li_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

  end do


!!!!!!!!!!!!!!!!!!
!!!  alpha<0   !!!
!!!!!!!!!!!!!!!!!!

 npts = 4

  alphamin=-0.1
  alphamax=-1.4*10._kp**(-3.)
  nalpha=40

!  w = 1._kp/3._kp
  w=0._kp

    do k=0,nalpha
      alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !logarithmic step
      alpha=-exp(log(-alphamin)/((log(-alphamin)/log(-alphamax))** &
            (real(k,kp)/real(nalpha,kp)))) !adapted step



  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = li_lnrhoreh_max(alpha,Pstar)

  print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax


  do i=0,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i,kp)/real(npts-1,kp)

       xstar = li_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

       eps1 = li_epsilon_one(xstar,alpha)
       eps2 = li_epsilon_two(xstar,alpha)
       eps3 = li_epsilon_three(xstar,alpha)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1*=',eps1

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('li_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('li_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

  end do


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('li_predic_summarized.dat') 
         nalpha=1000
         alphamin=0.002
         alphamax=10._kp**(5.)
         w=0._kp
         do j=1,nalpha
         alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))
         lnRhoReh = lnRhoNuc
         xstarA = li_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
         eps1A = li_epsilon_one(xstarA,alpha)
         eps2A = li_epsilon_two(xstarA,alpha)
         eps3A = li_epsilon_three(xstarA,alpha)
         nsA = 1._kp - 2._kp*eps1A - eps2A
         rA = 16._kp*eps1A
         lnRhoReh = li_lnrhoreh_max(alpha,Pstar)
         xstarB = li_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
         eps1B = li_epsilon_one(xstarB,alpha)
         eps2B = li_epsilon_two(xstarB,alpha)
         eps3B = li_epsilon_three(xstarB,alpha)
         nsB = 1._kp - 2._kp*eps1B - eps2B
         rB =16._kp*eps1B
         call livewrite('li_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
         enddo

end program limain
