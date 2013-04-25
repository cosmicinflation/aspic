!test the reheating derivation from slow-roll
program hf1imain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use hf1isr, only : hf1i_epsilon_one, hf1i_epsilon_two, hf1i_epsilon_three
  use hf1ireheat, only : hf1i_lnrhoreh_max, hf1i_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20
  integer ::nA1

  real(kp) :: A1,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer
  real(kp) ::A1min,A1max

  real(kp)  ::eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB


  nA1=10
  A1min=10._kp**(-3.)
  A1max=10._kp**(3.)

  Pstar = powerAmpScalar

!  w = 1._kp/3._kp
  w=0._kp

  call delete_file('hf1i_predic.dat')
  call delete_file('hf1i_nsr.dat')

  do j=0,nA1
    A1=A1min*(A1max/A1min)**(real(j,kp)/real(nA1,kp))


 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = hf1i_lnrhoreh_max(A1,Pstar)

  print *,'A1=',A1,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = hf1i_x_star(A1,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

       eps1 = hf1i_epsilon_one(xstar,A1)
       eps2 = hf1i_epsilon_two(xstar,A1)
       eps3 = hf1i_epsilon_three(xstar,A1)

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('hf1i_predic.dat',A1,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('hf1i_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('hf1i_predic_summarized.dat') 
         nA1=1000
         A1min=10._kp**(-3.)
         A1max=10._kp**(3.)
         w=0._kp
         do j=1,nA1
         A1=A1min*(A1max/A1min)**(real(j,kp)/real(nA1,kp))
         lnRhoReh = lnRhoNuc
         xstarA = hf1i_x_star(A1,w,lnRhoReh,Pstar,bfoldstar)
         eps1A = hf1i_epsilon_one(xstarA,A1)
         eps2A = hf1i_epsilon_two(xstarA,A1)
         eps3A = hf1i_epsilon_three(xstarA,A1)
         nsA = 1._kp - 2._kp*eps1A - eps2A
         rA = 16._kp*eps1A
         lnRhoReh = hf1i_lnrhoreh_max(A1,Pstar)
         xstarB = hf1i_x_star(A1,w,lnRhoReh,Pstar,bfoldstar)
         eps1B = hf1i_epsilon_one(xstarB,A1)
         eps2B = hf1i_epsilon_two(xstarB,A1)
         eps3B = hf1i_epsilon_three(xstarB,A1)
         nsB = 1._kp - 2._kp*eps1B - eps2B
         rB =16._kp*eps1B
         call livewrite('hf1i_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
         enddo



end program hf1imain
