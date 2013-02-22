!test the reheating derivation from slow-roll
program dwimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use dwisr, only : dwi_epsilon_one, dwi_epsilon_two, dwi_epsilon_three
  use dwireheat, only : dwi_lnrhoend, dwi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  integer :: Nphi0
 ! real(kp) :: phi0min=2._kp*sqrt(2._kp)
  real(kp) :: phi0min=7.6_kp
  real(kp) :: phi0max=10._kp**3

  real(kp) :: phi0,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB


  Nphi0=26

  Pstar = powerAmpScalar

  call delete_file('dwi_predic.dat')
  call delete_file('dwi_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp

 do j=0,Nphi0 

 phi0=phi0min+(phi0max-phi0min)*(real(j,kp)/real(Nphi0,kp)) !arithmetic step
 phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(Nphi0,kp)) !logarithmic step
 phi0=exp(log(phi0min)*(log(phi0max)/log(phi0min))**(real(j,kp)/real(Nphi0,kp))) !superlogarithmic step
 phi0=exp(exp(log(log(phi0min))*(log(log(phi0max))/log(log(phi0min)))**(real(j,kp)/real(Nphi0,kp)))) !ultralogarithmic step



  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = dwi_lnrhoend(phi0,Pstar)

  print *,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     

	xstar = dwi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)



       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar
 

       eps1 = dwi_epsilon_one(xstar,phi0)
       eps2 = dwi_epsilon_two(xstar,phi0)
       eps3 = dwi_epsilon_three(xstar,phi0)


       logErehGeV = log_energy_reheat_ingev(lnRhoReh)


       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('dwi_predic.dat',phi0,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('dwi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('dwi_predic_summarized.dat') 
         nphi0=1000
         phi0min=7.6_kp
         phi0max=10._kp**2
         w=0._kp
         do j=1,nphi0
         phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(nphi0,kp))
         lnRhoReh = lnRhoNuc
         xstarA = dwi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)
         eps1A = dwi_epsilon_one(xstarA,phi0)
         eps2A = dwi_epsilon_two(xstarA,phi0)
         eps3A = dwi_epsilon_three(xstarA,phi0)
         nsA = 1._kp - 2._kp*eps1A - eps2A
         rA = 16._kp*eps1A
         lnRhoReh = dwi_lnrhoend(phi0,Pstar)
         xstarB = dwi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)
         eps1B = dwi_epsilon_one(xstarB,phi0)
         eps2B = dwi_epsilon_two(xstarB,phi0)
         eps3B = dwi_epsilon_three(xstarB,phi0)
         nsB = 1._kp - 2._kp*eps1B - eps2B
         rB =16._kp*eps1B
         call livewrite('dwi_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
         enddo

 
end program dwimain
