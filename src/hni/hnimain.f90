!test the reheating derivation from slow-roll
program hnimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use hnisr, only : hni_epsilon_one, hni_epsilon_two,hni_epsilon_three
  use hnireheat, only : hni_lnrhoreh_max, hni_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use hnisr, only : hni_norm_potential, hni_x_endinf, hni_alphamin, hni_phizeromax
  use hnireheat, only : hni_x_rreh, hni_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh

  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nphi0

  real(kp) :: alpha,phi0,w,bfoldstar,alphamin,alphamax
  real(kp) :: phi0min,phi0max,phicutmin,phicutmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(:), allocatable ::phi0values

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!  allocate(phi0values(1:7))

!  phi0values(1)=5._kp
!  phi0values(2)=7._kp
!  phi0values(3)=10._kp
!  phi0values(4)=20._kp
!  phi0values(5)=50._kp
!  phi0values(6)=10._kp**(2.)

  npts = 15
  nalpha= 2
  nphi0 = 100

  phicutmin = 1._kp
  phicutmax = 200._kp
  
  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('hni_predic.dat')
  call delete_file('hni_nsr.dat')

  call aspicwrite_header('hni',labeps12,labnsr,labbfoldreh,(/'phi0   ','alpham1'/))



  alphamin = hni_alphamin(51._kp)
  alphamax = 1._kp-epsilon(1._kp)


  phi0min= phicutmin
  
  do k=0,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp))

     do j=0,nphi0

        phi0max = min(hni_phizeromax(alpha),phicutmax)
        
        phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(nphi0,kp))

           

        lnRhoRehMin = lnRhoNuc
        xEnd = hni_x_endinf(alpha,phi0)
        lnRhoRehMax = hni_lnrhoreh_max(alpha,phi0,xend,Pstar)

        print *,'alpha=',alpha,'phi0/Mp=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = hni_x_star(alpha,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)

           eps1 = hni_epsilon_one(xstar,alpha,phi0)
           eps2 = hni_epsilon_two(xstar,alpha,phi0)
           eps3 = hni_epsilon_three(xstar,alpha,phi0)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('hni_predic.dat',alpha,phi0,eps1,eps2,eps3,ns,r,Treh)

           call livewrite('hni_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/phi0,alpha-1._kp/))

        end do

     end do

  end do

  call aspicwrite_end()

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 0.999
  phi0 = 10.
  xEnd = hni_x_endinf(alpha,phi0)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = hni_x_rrad(alpha,phi0,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = hni_epsilon_one(xstar,alpha,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  hni_epsilon_one(xend,alpha,phi0)
     VendOverVstar = hni_norm_potential(xend,alpha,phi0)/hni_norm_potential(xstar,alpha,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = hni_x_rreh(alpha,phi0,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = hni_x_star(alpha,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo




end program hnimain
