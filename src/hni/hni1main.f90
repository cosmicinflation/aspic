!test the reheating derivation from slow-roll
program hni1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use hni1sr, only : hni_epsilon_one, hni_epsilon_two,hni_epsilon_three
  use hni1reheat, only : hni_lnrhoreh_max, hni_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use hni1sr, only : hni_norm_potential, hni_x_endinf, hni_alphamin, hni_phizeromax
  use hnir1eheat, only : hni_x_rreh, hni_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh

  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nphi0

  real(kp) :: alpha,f,w,bfoldstar,alphamin,alphamax
  real(kp) :: fmin,fmax,phicutmin,phicutmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(:), allocatable ::fvalues

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!  allocate(fvalues(1:7))

!  fvalues(1)=5._kp
!  fvalues(2)=7._kp
!  fvalues(3)=10._kp
!  fvalues(4)=20._kp
!  fvalues(5)=50._kp
!  fvalues(6)=10._kp**(2.)

  npts = 15
  nalpha= 2
  nf = 100

  phicutmin = 1._kp
  phicutmax = 200._kp
  
  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('hni1_predic.dat')
  call delete_file('hni1_nsr.dat')

  call aspicwrite_header('hni1',labeps12,labnsr,labbfoldreh,(/'f      ','alpham1'/))



  alphamin = hni_alphamin(51._kp)
  alphamax = 1._kp-epsilon(1._kp)


  fmin= phicutmin
  
  do k=0,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp))

     do j=0,nf

        fmax = min(hni_phizeromax(alpha),phicutmax)
        
        f=fmin*(fmax/fmin)**(real(j,kp)/real(nf,kp))

           

        lnRhoRehMin = lnRhoNuc
        xEnd = hni_x_endinf(alpha,f)
        lnRhoRehMax = hni_lnrhoreh_max(alpha,f,xend,Pstar)

        print *,'alpha=',alpha,'f/Mp=',f,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = hni_x_star(alpha,f,xend,w,lnRhoReh,Pstar,bfoldstar)

           eps1 = hni_epsilon_one(xstar,alpha,f)
           eps2 = hni_epsilon_two(xstar,alpha,f)
           eps3 = hni_epsilon_three(xstar,alpha,f)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('hni1_predic.dat',alpha,f,eps1,eps2,eps3,ns,r,Treh)

           call livewrite('hni1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/f,alpha-1._kp/))

        end do

     end do

  end do

  call aspicwrite_end()

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 0.999
  f = 10.
  xEnd = hni_x_endinf(alpha,f)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = hni_x_rrad(alpha,f,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = hni_epsilon_one(xstar,alpha,f)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  hni_epsilon_one(xend,alpha,f)
     VendOverVstar = hni_norm_potential(xend,alpha,f)/hni_norm_potential(xstar,alpha,f)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = hni_x_rreh(alpha,f,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = hni_x_star(alpha,f,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo




end program hni1main
