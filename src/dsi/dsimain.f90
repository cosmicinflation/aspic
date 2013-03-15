!test the reheating derivation from slow-roll
program dsimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use dsisr, only : dsi_epsilon_one, dsi_epsilon_two, dsi_epsilon_three,dsi_xinimin,dsi_xendmin,dsi_xendmax,dsi_mumax
  use dsireheat, only : dsi_lnrhoend, dsi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev


  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,l
  integer :: npts,nmu,nxEnd,np

  real(kp) :: p,mu,xEnd,w,bfoldstar,q
  real(kp) :: mumin,mumax,xEndmin,xEndmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp), dimension(:), allocatable ::pvalues,qvalues,muminvalues,mumaxvalues
  integer(kp), dimension(:), allocatable ::nmuvalues, nxEndvalues

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

  np=3
  allocate(pvalues(1:np))
  pvalues(1)=2._kp
  pvalues(2)=3._kp
  pvalues(3)=4._kp

  allocate(qvalues(1:np))
  qvalues(1)=8._kp
  qvalues(2)=8._kp
  qvalues(3)=8._kp

  allocate(nmuvalues(1:np))
  nmuvalues(1)=15
  nmuvalues(2)=15
  nmuvalues(3)=15

  allocate(mumaxvalues(1:np))
  mumaxvalues(1)=dsi_mumax(pvalues(1),qvalues(1),60._kp)
  mumaxvalues(2)=dsi_mumax(pvalues(2),qvalues(2),60._kp)
  mumaxvalues(3)=dsi_mumax(pvalues(3),qvalues(3),60._kp)

  allocate(muminvalues(1:np))
  muminvalues(1)=mumaxvalues(1)*10.**(-5.)
  muminvalues(2)=mumaxvalues(2)*10.**(-6.)
  muminvalues(3)=mumaxvalues(3)*10.**(-7.)

  allocate(nxEndvalues(1:np))
  nxEndvalues(1)=50
  nxEndvalues(2)=50
  nxEndvalues(3)=50

  do l=1,np
     p=pvalues(l)
     q=qvalues(l)
     nmu=nmuvalues(l)
     nxEnd=nxEndvalues(l)
     mumin=muminvalues(l)
     mumax=mumaxvalues(l)


  do j=0,nmu
     mu=mumin*(mumax/mumin)**(real(j,kp)/real(nmu,kp))

     xEndmin=dsi_xendmin(70._kp,p,mu)
     xEndmin=dsi_xendmin(60._kp,p,mu)
     xEndmax=dsi_xendmax(p,mu,q)
     print*,'xEndmin=',xEndmin,'xEndmax=',xEndmax

 
      do k=1,nxEnd
       xEnd=xEndmin*(xEndmax/xEndmin)**(real(k,kp)/real(nxEnd,kp)) 


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

end do




end program dsimain
