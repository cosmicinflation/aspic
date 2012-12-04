!test the reheating derivation from slow-roll
program cnbimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use cnbisr, only : cnbi_epsilon_one, cnbi_epsilon_two, cnbi_epsilon_three
  use cnbireheat, only : cnbi_lnrhoend, cnbi_x_star
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

  alphamin=10.**(-5)
  alphamax=0.2975_kp*0.99_kp

  Pstar = powerAmpScalar

  !w = 1._kp/3._kp
  w=0._kp

  call delete_file('cnbi_predic.dat')
  call delete_file('cnbi_nsr.dat')
  
  do j=1,nalpha
  alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))

  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = cnbi_lnrhoend(alpha,Pstar)

  print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = cnbi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)



       eps1 = cnbi_epsilon_one(xstar,alpha)
       eps2 = cnbi_epsilon_two(xstar,alpha)
       eps3 = cnbi_epsilon_three(xstar,alpha)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar, &
         'eps1star=',eps1,'eps2star=',eps2

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('cnbi_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('cnbi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do



end program cnbimain