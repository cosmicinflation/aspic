!test the reheating derivation from slow-roll
program bsusybimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use bsusybisr, only : bsusybi_epsilon_one, bsusybi_epsilon_two,bsusybi_epsilon_three,bsusybi_xendmax
  use bsusybireheat, only : bsusybi_lnrhoend, bsusybi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,ngamma,nxend

  real(kp) :: gammaBSUSYB,xend,w,bfoldstar,gammamin,gammamax,xendmin,xendmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer


  Pstar = powerAmpScalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          Calculates the prior space and               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ngamma=1000
  gammamin=10._kp**(-5._kp)
  gammamax=1._kp/sqrt(3._kp)

  call delete_file('bsusybi_xendmax.dat')
  do i=1,ngamma
       gammaBSUSYB=gammamin+(gammamax-gammamin)*(real(i,kp)/real(ngamma,kp))
       call livewrite('bsusybi_xendmax.dat',gammaBSUSYB,bsusybi_xendmax(gammaBSUSYB,-40._kp), &
       bsusybi_xendmax(gammaBSUSYB,-60._kp),bsusybi_xendmax(gammaBSUSYB,-80._kp))
  end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts = 20
  ngamma=10
  nxend=10

  gammamin=10._kp**(-3._kp)
  gammamax=1._kp/sqrt(3._kp)*0.5_kp

  w=0._kp
!  w = 1._kp/3._kp

  call delete_file('bsusybi_predic.dat')
  call delete_file('bsusybi_nsr.dat')


  do j=1,ngamma
    gammaBSUSYB=gammamin*(gammamax/gammamin)**(real(j,kp)/real(ngamma,kp))

  !Prior on xend

     xendmax=bsusybi_xendmax(gammaBSUSYB,-70._kp)
     xendmin=2._kp*xendmax


     do k=1,nxend
        xend=xendmin*(xendmax/xendmin)**(real(k,kp)/real(nxend,kp))
     

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = bsusybi_lnrhoend(gammaBSUSYB,xend,Pstar)

        print *,'gamma=',gammaBSUSYB,'xend=',xend,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = bsusybi_x_star(gammaBSUSYB,xend,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = bsusybi_epsilon_one(xstar,gammaBSUSYB)
           eps2 = bsusybi_epsilon_two(xstar,gammaBSUSYB)
           eps3 = bsusybi_epsilon_three(xstar,gammaBSUSYB)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('bsusybi_predic.dat',gammaBSUSYB,xend,xend/xendmax,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('bsusybi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
        end do

     end do

 end do




end program bsusybimain
