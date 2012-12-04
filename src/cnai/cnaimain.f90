!test the reheating derivation from slow-roll
program cnaimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use cnaisr, only : cnai_epsilon_one, cnai_epsilon_two, cnai_epsilon_three
  use cnaireheat, only : cnai_lnrhoend, cnai_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20
  integer :: nalpha = 300

  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::alphamin,alphamax

  alphamin=10.**(-4)
  alphamax=sqrt(0.5_kp*(sqrt(15._kp)-3._kp))

  Pstar = powerAmpScalar
  !w = 1._kp/3._kp
  w=0._kp

  call delete_file('cnai_predic.dat')
  call delete_file('cnai_nsr.dat')

  
  do j=1,nalpha
  alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))
 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = cnai_lnrhoend(alpha,Pstar)

  print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = cnai_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

       eps1 = cnai_epsilon_one(xstar,alpha)
       eps2 = cnai_epsilon_two(xstar,alpha)
       eps3 = cnai_epsilon_three(xstar,alpha)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar, &
         'eps1star=',eps1,'eps2star=',eps2

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('cnai_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('cnai_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do



end program cnaimain
