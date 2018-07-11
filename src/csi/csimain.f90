!test the reheating derivation from slow-roll
program csimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use csisr, only : csi_epsilon_one, csi_epsilon_two,csi_epsilon_three,csi_xendmax
  use csireheat, only : csi_lnrhoreh_max, csi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use csisr, only : csi_norm_potential
  use csireheat, only : csi_x_rreh, csi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nxend

  real(kp) :: alpha,xend,w,bfoldstar,alphamin,alphamax,xendmin,xendmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End

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

  npts = 8


  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('csi_predic.dat')
  call delete_file('csi_nsr.dat')

  call aspicwrite_header('csi',labeps12,labnsr,labbfoldreh,(/'xendomax','alpha   '/))
  
!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   alpha=0.001   !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

  alpha=10._kp**(-3.)

  !Prior on xend
  xendmax=csi_xendmax(65._kp,alpha)
  xendmin=-0.88*xendmax
  nxend=1000

  do k=1,nxend
     xend=xendmin+(xendmax-xendmin)*(real(k,kp)/real(nxend,kp))

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = csi_lnrhoreh_max(alpha,xend,Pstar)

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

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend/xendmax,alpha/))

     end do

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!     alpha=1     !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

  alpha=1._kp

  !Prior on xend
  xendmax=csi_xendmax(59._kp,alpha)
  xendmin=-1000._kp
  nxend=600

  do k=1,nxend
     xend=xendmin+(xendmax-xendmin)*(real(k,kp)/real(nxend,kp))

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = csi_lnrhoreh_max(alpha,xend,Pstar)

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

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend/xendmax,alpha/))

     end do

  end do


  call aspicwrite_end()

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 0.5
  xend = -20
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = csi_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = csi_epsilon_one(xstar,alpha)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  csi_epsilon_one(xend,alpha)
     VendOverVstar = csi_norm_potential(xend,alpha)/csi_norm_potential(xstar,alpha)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = csi_x_rreh(alpha,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = csi_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program csimain
