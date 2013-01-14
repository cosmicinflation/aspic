!test the reheating derivation from slow-roll
program cndimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use cndisr, only : cndi_epsilon_one, cndi_epsilon_two, cndi_epsilon_three, cndi_xend_max
  use cndireheat, only : cndi_lnrhoend, cndi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev


  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nbeta,nxend

  real(kp) :: alpha,beta,xend,w,bfoldstar,alphamin,alphamax,betamin,betamax,xendmin,xendmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(10) :: alphavalues
  integer, dimension(10) :: nxendvalues


  Pstar = powerAmpScalar



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    w=0._kp
!  w = 1._kp/3._kp

  call delete_file('cndi_predic.dat')
  call delete_file('cndi_nsr.dat')



  npts = 20

  beta=5._kp
  beta=0.1_kp

  if (beta .eq. 0.1_kp) then
  nalpha=6
  alphavalues(1)=0.1_kp
  alphavalues(2)=0.2_kp
  alphavalues(3)=0.3_kp
  alphavalues(4)=0.4_kp
  alphavalues(5)=0.5_kp
  alphavalues(6)=0.6_kp

  nxendvalues(1)=1100
  nxendvalues(2)=800
  nxendvalues(3)=400
  nxendvalues(4)=100
  nxendvalues(5)=150
  nxendvalues(6)=100
  endif

  if (beta .eq. 5._kp) then
  nalpha=4
  alphavalues(1)=0.01_kp
  alphavalues(2)=0.07_kp
  alphavalues(3)=0.1_kp
  alphavalues(4)=0.12_kp

  nxendvalues(1)=1000
  nxendvalues(2)=200
  nxendvalues(3)=200
  nxendvalues(4)=200

  endif

  do j=1,nalpha
       alpha=alphavalues(j)
       nxend=nxendvalues(j)

  xendmin = sqrt(epsilon(1._kp)/2._kp)*(1._kp+beta)/(alpha**2*beta) !to avoid  epsilon_1 < numaccuracy errors
  xendmax = cndi_xend_max(70._kp,alpha,beta)
  xendmin = xendmax/1000._kp


  do k=0,nxend
       xend=xendmin*(xendmax/xendmin)**(real(k,kp)/real(nxend,kp)) !log step
     

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = cndi_lnrhoend(alpha,beta,xend,Pstar)

  print *,'alpha=',alpha,'beta=',beta,'xend=',xend,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax


  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = cndi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)

       eps1 = cndi_epsilon_one(xstar,alpha,beta)
       eps2 = cndi_epsilon_two(xstar,alpha,beta)
       eps3 = cndi_epsilon_three(xstar,alpha,beta)


       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('cndi_predic.dat',alpha,beta,xend,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('cndi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do
  end do
end do




end program cndimain
