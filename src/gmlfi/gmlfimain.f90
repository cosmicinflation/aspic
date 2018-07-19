!test the reheating derivation from slow-roll
program gmlfimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use gmlfisr, only : gmlfi_epsilon_one, gmlfi_epsilon_two, gmlfi_epsilon_three
  use gmlfireheat, only : gmlfi_lnrhoreh_max, gmlfi_x_star
  use infinout, only : delete_file, livewrite, has_not_shifted
  use srreheat, only : log_energy_reheat_ingev

  use gmlfisr, only : gmlfi_norm_potential, gmlfi_x_endinf
  use gmlfireheat, only : gmlfi_x_rreh, gmlfi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,np,nq

  real(kp) :: alpha,p,q,w,bfoldstar,alphamin,alphamax,pmin,pmax,qmin,qmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::x,xmin,xmax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('gmlfi_predic.dat')
  call delete_file('gmlfi_nsr.dat')

  call aspicwrite_header('gmlfi',labeps12,labnsr,labbfoldreh,(/'alpha','q    ','p    '/))
  
  w=0._kp
  !  w = 1._kp/3._kp

  npts = 10

!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    p=2  &  q=1  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

  p=2.
  q=1.
  alphamin=10._kp**(-3._kp)
  alphamax=10._kp**(3._kp)
  nalpha=100

  do j=0,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))

     lnRhoRehMin = lnRhoNuc
     xEnd = gmlfi_x_endinf(p,q,alpha)
     lnRhoRehMax = gmlfi_lnrhoreh_max(p,q,alpha,xend,Pstar)


     print *,'alpha,p,q=',alpha,p,q,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = gmlfi_x_star(p,q,alpha,xend,w,lnRhoReh,Pstar,bfoldstar)


        eps1 = gmlfi_epsilon_one(xstar,p,q,alpha)
        eps2 = gmlfi_epsilon_two(xstar,p,q,alpha)
        eps3 = gmlfi_epsilon_three(xstar,p,q,alpha)


        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1   

        call livewrite('gmlfi_predic.dat',p,q,alpha,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('gmlfi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,q,p/))
        
     end do


  end do


!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    p=2  &  q=3  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!


  p=2.
  q=3.
  alphamin=10._kp**(-8._kp)
  alphamax=10._kp**(1._kp)
  nalpha=100

  do j=0,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))

     lnRhoRehMin = lnRhoNuc
     xEnd = gmlfi_x_endinf(p,q,alpha)
     lnRhoRehMax = gmlfi_lnrhoreh_max(p,q,alpha,xend,Pstar)


     print *,'alpha,p,q=',alpha,p,q,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = gmlfi_x_star(p,q,alpha,xend,w,lnRhoReh,Pstar,bfoldstar)


        eps1 = gmlfi_epsilon_one(xstar,p,q,alpha)
        eps2 = gmlfi_epsilon_two(xstar,p,q,alpha)
        eps3 = gmlfi_epsilon_three(xstar,p,q,alpha)


        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1


        call livewrite('gmlfi_predic.dat',p,q,alpha,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('gmlfi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,q,p/))

     end do


  end do


!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    p=3  &  q=2  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!

  p=3.
  q=2.
  alphamin=10._kp**(-6._kp)
  alphamax=10._kp**(3._kp)
  nalpha=100


  do j=0,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))

     lnRhoRehMin = lnRhoNuc
     xEnd = gmlfi_x_endinf(p,q,alpha)
     lnRhoRehMax = gmlfi_lnrhoreh_max(p,q,alpha,xend,Pstar)


     print *,'alpha,p,q=',alpha,p,q,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = gmlfi_x_star(p,q,alpha,xend,w,lnRhoReh,Pstar,bfoldstar)


        eps1 = gmlfi_epsilon_one(xstar,p,q,alpha)
        eps2 = gmlfi_epsilon_two(xstar,p,q,alpha)
        eps3 = gmlfi_epsilon_three(xstar,p,q,alpha)


        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1



        call livewrite('gmlfi_predic.dat',p,q,alpha,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('gmlfi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha,q,p/))

     end do


  end do

  call aspicwrite_end()
  
  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  p=2
  q=3
  alpha = 0.01
  xEnd = gmlfi_x_endinf(p,q,alpha)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = gmlfi_x_rrad(p,q,alpha,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = gmlfi_epsilon_one(xstar,p,q,alpha)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = gmlfi_x_endinf(p,q,alpha)
     eps1end =  gmlfi_epsilon_one(xend,p,q,alpha)
     VendOverVstar = gmlfi_norm_potential(xend,p,q,alpha)/gmlfi_norm_potential(xstar,p,q,alpha)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = gmlfi_x_rreh(p,q,alpha,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = gmlfi_x_star(p,q,alpha,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program gmlfimain
