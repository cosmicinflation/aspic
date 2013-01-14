!test the reheating derivation from slow-roll
program cncimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use cncisr, only : cnci_epsilon_one, cnci_epsilon_two,cnci_epsilon_three,cnci_xendmin,cnci_x_epsOne_equals_one
  use cncireheat, only : cnci_lnrhoend, cnci_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nxend

  real(kp) :: alpha,xend,w,bfoldstar,alphamin,alphamax,xendmin,xendmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r
  real(kp), dimension(10) ::alphavalues

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

  call delete_file('cnci_xendmax.dat')
  do i=1,nalpha
       alpha=alphamin+(alphamax-alphamin)*(real(i,kp)/real(nalpha,kp))

       call livewrite('cnci_xendmin.dat',alpha,cnci_xendmin(alpha,-40._kp), &
       cnci_xendmin(alpha,-60._kp),cnci_xendmin(alpha,-80._kp))
  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts = 20
  nxend=20


  w=0._kp
!  w = 1._kp/3._kp

  call delete_file('cnci_predic.dat')
  call delete_file('cnci_nsr.dat')

  alphavalues(3)=10._kp**(-3.)
  alphavalues(2)=10._kp**(-1.)
  alphavalues(1)=0.2_kp

  do j=1,3
   alpha=alphavalues(j)



  !Prior on xend
     if (alpha .eq. 10._kp**(-3.)) then
     xendmin=cnci_xendmin(alpha,-55._kp)
     xendmax=30._kp*xendmin
     nxend=500
     endif
     if (alpha .eq. 10._kp**(-1.)) then
     xendmin=cnci_xendmin(alpha,-58._kp)
     xendmax=10._kp*xendmin
     nxend=500
     endif
     if (alpha .eq. 0.2_kp) then
     xendmin=cnci_xendmin(alpha,-60._kp)
     xendmax=2._kp*xendmin
     nxend=500
     endif



     do k=1,nxend
        xend=xendmin+(xendmax-xendmin)*(real(k,kp)/real(nxend,kp))

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = cnci_lnrhoend(alpha,xend,Pstar)

        print *,'alpha=',alpha,'xend=',xend,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = cnci_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = cnci_epsilon_one(xstar,alpha)
           eps2 = cnci_epsilon_two(xstar,alpha)
           eps3 = cnci_epsilon_three(xstar,alpha)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('cnci_predic.dat',alpha,xend,xendmin,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('cnci_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
        end do

     end do

 end do




end program cncimain
