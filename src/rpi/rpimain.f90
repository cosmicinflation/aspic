!test the reheating derivation from slow-roll
program rpimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rpisr, only : rpi_epsilon_one, rpi_epsilon_two, rpi_epsilon_three
  use rpireheat, only : rpi_lnrhoend, rpi_y_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  integer :: Np=1
  real(kp) :: pmin=1._kp
  real(kp) :: pmax=10._kp

  real(kp) :: p,w,bfoldstar
  real(kp) :: lnRhoReh,ystar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer



  Pstar = powerAmpScalar

  call delete_file('rpi_predic.dat')
  call delete_file('rpi_nsr.dat')


!  w = 1._kp/3._kp
  w=0._kp

 do j=0,Np 
 p=pmin+(pmax-pmin)*(real(j,kp)/Np)
 w=(1._kp-p)/(3._kp-1._kp)


  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = rpi_lnrhoend(p,Pstar)

  print *,'p=',p,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     

	ystar = rpi_y_star(p,w,lnRhoReh,Pstar,bfoldstar)



       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'ystar=',ystar
 

       eps1 = rpi_epsilon_one(ystar,p)
       eps2 = rpi_epsilon_two(ystar,p)
       eps3 = rpi_epsilon_three(ystar,p)


       logErehGeV = log_energy_reheat_ingev(lnRhoReh)


       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('rpi_predic.dat',p,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('rpi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

 end do

 

end program rpimain
