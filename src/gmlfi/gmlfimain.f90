!test the reheating derivation from slow-roll
program gmlfimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use gmlfisr, only : gmlfi_epsilon_one, gmlfi_epsilon_two, gmlfi_epsilon_three
  use gmlfireheat, only : gmlfi_lnrhoend, gmlfi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev


  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,np,nq

  real(kp) :: alpha,p,q,w,bfoldstar,alphamin,alphamax,pmin,pmax,qmin,qmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::x,xmin,xmax


  Pstar = powerAmpScalar



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('gmlfi_predic.dat')
  call delete_file('gmlfi_nsr.dat')

  w=0._kp
!  w = 1._kp/3._kp

  npts = 20

  p=2.7
  q=4.8
  
  nalpha=20

  alphamin=10._kp**(-3._kp)
  alphamax=10._kp**(3._kp)

  do j=0,nalpha
    alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))
 
    lnRhoRehMin = lnRhoNuc
    lnRhoRehMax = gmlfi_lnrhoend(p,q,alpha,Pstar)


    print *,'alpha,p,q=',alpha,p,q,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = gmlfi_x_star(p,q,alpha,w,lnRhoReh,Pstar,bfoldstar)


       eps1 = gmlfi_epsilon_one(xstar,p,q,alpha)
       eps2 = gmlfi_epsilon_two(xstar,p,q,alpha)
       eps3 = gmlfi_epsilon_three(xstar,p,q,alpha)


       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('gmlfi_predic.dat',alpha,p,q,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('gmlfi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do


 end do



end program gmlfimain
