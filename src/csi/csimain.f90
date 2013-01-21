!test the reheating derivation from slow-roll
program csimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use csisr, only : csi_epsilon_one, csi_epsilon_two,csi_epsilon_three,csi_xendmax
  use csireheat, only : csi_lnrhoend, csi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nxend

  real(kp) :: alpha,xend,w,bfoldstar,alphamin,alphamax,xendmin,xendmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer


  Pstar = powerAmpScalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!          Calculates the prior space and               !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  nalpha=1000
  alphamin=10._kp**(-5._kp)
  alphamax=1._kp

  call delete_file('csi_xendmax.dat')
  do i=1,nalpha
       alpha=alphamin+(alphamax-alphamin)*(real(i,kp)/real(nalpha,kp))

       call livewrite('csi_xendmax.dat',alpha,csi_xendmax(40._kp,alpha), &
       csi_xendmax(60._kp,alpha),csi_xendmax(80._kp,alpha))
  end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts = 20
!  nalpha=10
  nxend=20

!  alphamin=10._kp**(-3._kp)
!  alphamax=10._kp**(3._kp)

  w=0._kp
!  w = 1._kp/3._kp

  call delete_file('csi_predic.dat')
  call delete_file('csi_nsr.dat')


!  do j=1,nalpha
!    alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))

  alpha=10._kp**(-3.)
  alpha=1._kp

  !Prior on xend
     if (alpha .eq. 10._kp**(-3.)) then
     xendmax=csi_xendmax(55._kp,alpha)
     xendmin=-xendmax
     nxend=500
     endif
     if (alpha .eq. 1._kp) then
     xendmax=csi_xendmax(55._kp,alpha)
     xendmin=-1000._kp
     nxend=400
     endif


     do k=1,nxend
        xend=xendmin+(xendmax-xendmin)*(real(k,kp)/real(nxend,kp))

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = csi_lnrhoend(alpha,xend,Pstar)

        print *,'alpha=',alpha,'xend=',xend,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = csi_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = csi_epsilon_one(xstar,alpha)
           eps2 = csi_epsilon_two(xstar,alpha)
           eps3 = csi_epsilon_three(xstar,alpha)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('csi_predic.dat',alpha,xend,abs(1._kp-xend/xendmax),eps1,eps2,eps3,r,ns,Treh)

           call livewrite('csi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
        end do

!     end do

 end do




end program csimain
