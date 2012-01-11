!test the reheating derivation from slow-roll
program cwimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use cwisr, only : cwi_epsilon_one, cwi_epsilon_two, cwi_epsilon_three
  use cwireheat, only : cwi_lnrhoend, cwi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  integer :: Nalpha=10
  real(kp) :: alphamin=0.00000001
  real(kp) :: alphamax=0.01

  real(kp) :: alpha,Q,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer



  Pstar = powerAmpScalar

  call delete_file('cwi_predic.dat')
  call delete_file('cwi_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp

 do j=0,Nalpha 
 alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/Nalpha)
 Q=(4._kp*exp(1._kp)/alpha)**(1._kp/4._kp)


  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = cwi_lnrhoend(alpha,Q,Pstar)

  print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     

	xstar = cwi_x_star(alpha,Q,w,lnRhoReh,Pstar,bfoldstar)



       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar
 

       eps1 = cwi_epsilon_one(xstar,alpha,Q)
       eps2 = cwi_epsilon_two(xstar,alpha,Q)
       eps3 = cwi_epsilon_three(xstar,alpha,Q)


       logErehGeV = log_energy_reheat_ingev(lnRhoReh)


       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('cwi_predic.dat',alpha,Q,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('cwi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

 

end program cwimain
