!test the reheating derivation from slow-roll
program esimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use esisr, only : esi_epsilon_one, esi_epsilon_two, esi_epsilon_three
  use esireheat, only : esi_lnrhoend, esi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i
  integer :: npts = 20

  real(kp) :: q,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  Pstar = powerAmpScalar

!  call delete_file('esi_predic.dat')
!  call delete_file('esi_nsr.dat')

  q = 3.5_kp*10._kp**(0.) 
!  w = 1._kp/3._kp
  w=0._kp
 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = esi_lnrhoend(q,Pstar)

  print *,'q=',q,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = esi_x_star(q,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

       eps1 = esi_epsilon_one(xstar,q)
       eps2 = esi_epsilon_two(xstar,q)
       eps3 = esi_epsilon_three(xstar,q)


       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('esi_predic.dat',q,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('esi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do



end program esimain
