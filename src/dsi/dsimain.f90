!test the reheating derivation from slow-roll
program dsimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use dsisr, only : dsi_epsilon_one, dsi_epsilon_two, dsi_epsilon_three,dsi_xinimin
  use dsireheat, only : dsi_lnrhoend, dsi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev


  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nmu,nxEnd

  real(kp) :: p,mu,xEnd,w,bfoldstar
  real(kp) :: mumin,mumax,xEndmin,xEndmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer


  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('dsi_predic.dat')
  call delete_file('dsi_nsr.dat')

  npts = 20

  w=0._kp
!  w = 1._kp/3._kp

  p=1.5_kp

  nmu=20
  nxEnd=10

  mumin=10**(-1._kp)
  mumax=10**(1._kp)

  do j=0,nmu
     mu=mumin*(mumax/mumin)**(real(j,kp)/real(nmu,kp))

!     mu=0.1  !!!!!!!!!!!!!!

     xEndmin=dsi_xinimin(p,mu)*10.
     xEndmax=100._kp*xEndmin

     do k=1,nxEnd
        xEnd=xEndmin*(xEndmax/xEndmin)**(real(k,kp)/real(nxEnd,kp))

!        xEnd=10.*(p/(sqrt(2._kp)*mu))**(1._kp/(p+1._kp))  !!!!!!!!!!!!!!!!!!!!!!!


        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = dsi_lnrhoend(p,mu,xEnd,Pstar)

        print *,'p=',p,'mu=',mu,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = dsi_x_star(p,mu,xEnd,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = dsi_epsilon_one(xstar,p,mu)
           eps2 = dsi_epsilon_two(xstar,p,mu)
           eps3 = dsi_epsilon_three(xstar,p,mu)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('dsi_predic.dat',p,mu,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('dsi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do




end program dsimain
