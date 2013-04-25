!test the reheating derivation from slow-roll
program mhimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use mhisr, only : mhi_epsilon_one, mhi_epsilon_two, mhi_epsilon_three,mhi_x_endinf
  use mhireheat, only : mhi_lnrhoreh_max, mhi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh, mu, mumin, mumax

  integer :: i,j
  integer :: npts = 20, nj=20

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  integer :: nmu



  Pstar = powerAmpScalar

  call delete_file('mhi_predic.dat')
  call delete_file('mhi_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp

  mumin=(10._kp)**(-1)
  mumax=(10._kp)**2
  
  do j=0,nj
    mu=mumin*(mumax/mumin)**(real(j,kp)/real(nj,kp))


  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = mhi_lnrhoreh_max(mu,Pstar)

  print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  print*,'mu=',mu,'xEnd=',mhi_x_endinf(mu)

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     

	xstar = mhi_x_star(mu,w,lnRhoReh,Pstar,bfoldstar)



       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',mhi_epsilon_one(xstar,mu)
 

       eps1 = mhi_epsilon_one(xstar,mu)
       eps2 = mhi_epsilon_two(xstar,mu)
       eps3 = mhi_epsilon_three(xstar,mu)


       logErehGeV = log_energy_reheat_ingev(lnRhoReh)


       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('mhi_predic.dat',mu,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('mhi_nsr.dat',mu,ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('mhi_predic_summarized.dat') 
         nmu=1000
         mumin=(10._kp)**(-1)
         mumax=(10._kp)**2
         w=0._kp
         do j=1,nmu
         mu=mumin*(mumax/mumin)**(real(j,kp)/real(nmu,kp))
         lnRhoReh = lnRhoNuc
         xstarA = mhi_x_star(mu,w,lnRhoReh,Pstar,bfoldstar)
         eps1A = mhi_epsilon_one(xstarA,mu)
         eps2A = mhi_epsilon_two(xstarA,mu)
         eps3A = mhi_epsilon_three(xstarA,mu)
         nsA = 1._kp - 2._kp*eps1A - eps2A
         rA = 16._kp*eps1A
         lnRhoReh = mhi_lnrhoreh_max(mu,Pstar)
         xstarB = mhi_x_star(mu,w,lnRhoReh,Pstar,bfoldstar)
         eps1B = mhi_epsilon_one(xstarB,mu)
         eps2B = mhi_epsilon_two(xstarB,mu)
         eps3B = mhi_epsilon_three(xstarB,mu)
         nsB = 1._kp - 2._kp*eps1B - eps2B
         rB =16._kp*eps1B
         call livewrite('mhi_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
         enddo

end program mhimain
