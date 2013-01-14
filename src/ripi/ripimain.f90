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
  integer :: nalpha

  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::alphamin,alphamax

  real(kp) :: eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB


  nalpha = 100

  alphamin=10.**(-4.)
  alphamax=10.**(0.)

  Pstar = powerAmpScalar

  w=0._kp
  !w = 1._kp/3._kp

  call delete_file('ripi_predic.dat')
  call delete_file('ripi_nsr.dat')


  do j=1,nalpha

  alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))


  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = ripi_lnrhoend(alpha,Pstar)

  print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = ripi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

       eps1 = ripi_epsilon_one(xstar,alpha)
       eps2 = ripi_epsilon_two(xstar,alpha)
       eps3 = ripi_epsilon_three(xstar,alpha)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar, &
         'eps1star=',eps1,'eps2star=',eps2

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('ripi_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('ripi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('ripi_predic_summarized.dat') 
         nalpha=1000
         alphamin=10.**(-4.)
         alphamax=10.**(0.)
         w=0._kp
         do j=1,nalpha
         alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))
         lnRhoReh = lnRhoNuc
         xstarA = ripi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
         eps1A = ripi_epsilon_one(xstarA,alpha)
         eps2A = ripi_epsilon_two(xstarA,alpha)
         eps3A = ripi_epsilon_three(xstarA,alpha)
         nsA = 1._kp - 2._kp*eps1A - eps2A
         rA = 16._kp*eps1A
         lnRhoReh = ripi_lnrhoend(alpha,Pstar)
         xstarB = ripi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
         eps1B = ripi_epsilon_one(xstarB,alpha)
         eps2B = ripi_epsilon_two(xstarB,alpha)
         eps3B = ripi_epsilon_three(xstarB,alpha)
         nsB = 1._kp - 2._kp*eps1B - eps2B
         rB =16._kp*eps1B
         call livewrite('ripi_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
         enddo



end program ripimain
