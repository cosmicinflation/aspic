!test the reheating derivation from slow-roll
program rcmimain
  use infprec, only : kp, pi
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rcmisr, only : rcmi_norm_potential
  use rcmisr, only : rcmi_epsilon_one, rcmi_epsilon_two, rcmi_epsilon_three
  use rcmireheat, only : rcmi_lnrhoreh_max 
  use rcmireheat, only : rcmi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use rcmisr, only : rcmi_norm_potential, rcmi_x_endinf
  use rcmireheat, only : rcmi_x_rreh, rcmi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  real(kp) :: alpha,w,bfoldstar,Vstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(1:7) ::alphavalues

  real(kp)  ::alphamin,alphamax,eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  integer :: nalpha
  
  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend
  

  alphavalues(1)=(10._kp)**(-6.)
  alphavalues(2)=5.*(10._kp)**(-5.)
  alphavalues(3)=(10._kp)**(-4.)
  alphavalues(4)=2.*(10._kp)**(-4.)
  alphavalues(5)=3.5*(10._kp)**(-4.)
  alphavalues(6)=4.5*(10._kp)**(-4.)
  alphavalues(7)=5.7*(10._kp)**(-4.)

  Pstar = powerAmpScalar
  w = 0._kp

  call delete_file('rcmi_predic.dat')
  call delete_file('rcmi_nsr.dat')

  !  alphamin=alphavalues(1)
  !  alphamax=alphavalues(7)
  !  do j=1,40
  !  alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/40.)


  do j=1,size(alphavalues)  
     alpha=alphavalues(j)


     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = rcmi_lnrhoreh_max(alpha,Pstar)

     print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = rcmi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,'bfoldstar= ',bfoldstar, 'xstar=',xstar

        eps1 = rcmi_epsilon_one(xstar,alpha)
        eps2 = rcmi_epsilon_two(xstar,alpha)
        eps3 = rcmi_epsilon_three(xstar,alpha)

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh =  10._kp**(logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp))

        Vstar=(3._kp-eps1)*(Pstar*8*pi*pi*eps1)
        print *,'M4=',Vstar/rcmi_norm_potential(xstar,alpha)

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('rcmi_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('rcmi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('rcmi_predic_summarized.dat') 
  nalpha=1000
  alphamin=10._kp**(-6.)
  alphamax=10._kp**(-3.)
  do j=1,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))
     lnRhoReh = lnRhoNuc
     xstarA = rcmi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     eps1A = rcmi_epsilon_one(xstarA,alpha)
     eps2A = rcmi_epsilon_two(xstarA,alpha)
     eps3A = rcmi_epsilon_three(xstarA,alpha)
     nsA = 1._kp - 2._kp*eps1A - eps2A
     rA = 16._kp*eps1A
     lnRhoReh = rcmi_lnrhoreh_max(alpha,Pstar)
     xstarB = rcmi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     eps1B = rcmi_epsilon_one(xstarB,alpha)
     eps2B = rcmi_epsilon_two(xstarB,alpha)
     eps3B = rcmi_epsilon_three(xstarB,alpha)
     nsB = 1._kp - 2._kp*eps1B - eps2B
     rB =16._kp*eps1B
     call livewrite('rcmi_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
  enddo


write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 5e-4
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rcmi_x_rrad(alpha,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rcmi_epsilon_one(xstar,alpha)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = rcmi_x_endinf(alpha)
     eps1end =  rcmi_epsilon_one(xend,alpha)
     VendOverVstar = rcmi_norm_potential(xend,alpha)/rcmi_norm_potential(xstar,alpha)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rcmi_x_rreh(alpha,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rcmi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program rcmimain
