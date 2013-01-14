!test the reheating derivation from slow-roll
program rgimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rgisr, only : rgi_epsilon_one, rgi_epsilon_two, rgi_epsilon_three
  use rgireheat, only : rgi_lnrhoend, rgi_lnrhoreh_fromepsilon 
  use rgireheat, only : rgi_xp_fromepsilon, rgi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  integer :: Nalpha
  real(kp) :: alphamin=0.00001
  real(kp) :: alphamax=10000._kp

  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer


  real(kp) :: eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB

  Nalpha = 20

  Pstar = powerAmpScalar

  call delete_file('rgi_predic.dat')
  call delete_file('rgi_nsr.dat')

 w=0.


 do j=0,Nalpha
    
     alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(Nalpha,kp))
 
    lnRhoRehMin = lnRhoNuc
    lnRhoRehMax = rgi_lnrhoend(alpha,Pstar)

    print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

    do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = rgi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

       eps1 = rgi_epsilon_one(xstar,alpha)
       eps2 = rgi_epsilon_two(xstar,alpha)
       eps3 = rgi_epsilon_three(xstar,alpha)

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('rgi_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('rgi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('rgi_predic_summarized.dat') 
         nalpha=1000
         alphamin=10._kp**(-5.)
         alphamax=10._kp**(4.)
         w=0._kp
         do j=1,nalpha
         alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))
         lnRhoReh = lnRhoNuc
         xstarA = rgi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
         eps1A = rgi_epsilon_one(xstarA,alpha)
         eps2A = rgi_epsilon_two(xstarA,alpha)
         eps3A = rgi_epsilon_three(xstarA,alpha)
         nsA = 1._kp - 2._kp*eps1A - eps2A
         rA = 16._kp*eps1A
         lnRhoReh = rgi_lnrhoend(alpha,Pstar)
         xstarB = rgi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
         eps1B = rgi_epsilon_one(xstarB,alpha)
         eps2B = rgi_epsilon_two(xstarB,alpha)
         eps3B = rgi_epsilon_three(xstarB,alpha)
         nsB = 1._kp - 2._kp*eps1B - eps2B
         rB =16._kp*eps1B
         call livewrite('rgi_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
         enddo


end program rgimain
