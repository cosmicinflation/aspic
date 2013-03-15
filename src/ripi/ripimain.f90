!test the reheating derivation from slow-roll
program ripimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ripisr, only : ripi_epsilon_one, ripi_epsilon_two, ripi_epsilon_three
  use ripireheat, only : ripi_lnrhoend, ripi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20
  integer :: nphi0

  real(kp) :: phi0,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::phi0min,phi0max

  real(kp) :: eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB


  nphi0 = 50

  phi0min=10.**(-5.)
  phi0max=10.**(3.)

  Pstar = powerAmpScalar

  w=0._kp
  !w = 1._kp/3._kp

  call delete_file('ripi_predic.dat')
  call delete_file('ripi_nsr.dat')


  do j=1,nphi0

  phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(nphi0,kp))


  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = ripi_lnrhoend(phi0,Pstar)

  print *,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = ripi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)

       eps1 = ripi_epsilon_one(xstar,phi0)
       eps2 = ripi_epsilon_two(xstar,phi0)
       eps3 = ripi_epsilon_three(xstar,phi0)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar, &
         'eps1star=',eps1,'eps2star=',eps2

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('ripi_predic.dat',phi0,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('ripi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('ripi_predic_summarized.dat') 
         nphi0=1000
         phi0min=10.**(-5.)
         phi0max=10.**(-3.)
         w=0._kp
         do j=1,nphi0
         phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(nphi0,kp))
         lnRhoReh = lnRhoNuc
         xstarA = ripi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)
         eps1A = ripi_epsilon_one(xstarA,phi0)
         eps2A = ripi_epsilon_two(xstarA,phi0)
         eps3A = ripi_epsilon_three(xstarA,phi0)
         nsA = 1._kp - 2._kp*eps1A - eps2A
         rA = 16._kp*eps1A
         lnRhoReh = ripi_lnrhoend(phi0,Pstar)
         xstarB = ripi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)
         eps1B = ripi_epsilon_one(xstarB,phi0)
         eps2B = ripi_epsilon_two(xstarB,phi0)
         eps3B = ripi_epsilon_three(xstarB,phi0)
         nsB = 1._kp - 2._kp*eps1B - eps2B
         rB =16._kp*eps1B
         call livewrite('ripi_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
         enddo



end program ripimain
