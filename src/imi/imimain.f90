!test the reheating derivation from slow-roll
program imimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use imisr, only : imi_epsilon_one, imi_epsilon_two, imi_epsilon_three, imi_xendmin
  use imireheat, only : imi_lnrhoend, imi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,nxend
  integer :: npts = 20

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: pmin,pmax,p,xendmin,xendmax,xend

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  Pstar = powerAmpScalar

  nxend=30

  w=0._kp

  call delete_file('imi_predic.dat')
  call delete_file('imi_nsr.dat')

  pmin=1.
  pmax=6.

  do j=0,int(pmax-pmin)

  p = pmin+real(j,kp)

  xendmin=imi_xendmin(65._kp,p)
  xendmax=100._kp*xendmin

  do k=0,nxend
    
    xend=xendmin*(xendmax/xendmin)**(real(k,kp)/real(nxend,kp))
 
    lnRhoRehMin = lnRhoNuc
    lnRhoRehMax = imi_lnrhoend(p,xend,Pstar)

    print *,'p=',p,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

    do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = imi_x_star(p,xend,w,lnRhoReh,Pstar,bfoldstar)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar

       eps1 = imi_epsilon_one(xstar,p)
       eps2 = imi_epsilon_two(xstar,p)
       eps3 = imi_epsilon_three(xstar,p)

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('imi_predic.dat',p,xend,eps1,eps2,eps3,r,ns,Treh)

  
    end do

  end do

 end do

end program imimain
