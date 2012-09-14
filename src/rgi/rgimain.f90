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

  integer :: Nalpha = 20
  real(kp) :: alphamin=0.00001
  real(kp) :: alphamax=10000._kp

  real(kp) :: alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

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

end program rgimain
