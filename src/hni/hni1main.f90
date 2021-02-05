!test the reheating derivation from slow-roll
program hni1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar


  use hni1sr, only : hni1_norm_potential, hni1_x_endinf
  use hni1sr, only : hni1_epsilon_one, hni1_epsilon_two,hni1_epsilon_three
  use hni1sr, only : hni1_alphamin, hni1_fmax, hni1_numacc_efoldmax
  
  use hni1reheat, only : hni1_x_rreh, hni1_x_rrad
  use hni1reheat, only : hni1_lnrhoreh_max, hni1_x_star

  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat
  use srreheat, only : log_energy_reheat_ingev

  use infinout, only : delete_file, livewrite
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh

  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nf

  real(kp) :: alpha,f,w,bfoldstar,alphamin,alphamax
  real(kp) :: fmin,fmax,fcutmin,fcutmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(:), allocatable ::fvalues

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: efold, efoldmax
  
  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




  npts = 15
  nalpha= 4
  nf = 100

  efold = 120._kp
  
  fcutmin = 0.1_kp
  fcutmax = 100._kp
  
  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('hni1_predic.dat')
  call delete_file('hni1_nsr.dat')

  call aspicwrite_header('hni1',labeps12,labnsr,labbfoldreh,(/'f      ','1malpha'/))



  alphamin = hni1_alphamin(51._kp)
  alphamax = 0.99


  fmin= fcutmin
  
  do k=1,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(k-1,kp)/real(nalpha-1,kp))

     do j=1,nf

        fmax = min(0.99*hni1_fmax(alpha),fcutmax)
        
        f=fmin*(fmax/fmin)**(real(j-1,kp)/real(nf-1,kp))


        lnRhoRehMin = lnRhoNuc
        xEnd = hni1_x_endinf(alpha,f)

        efoldmax = hni1_numacc_efoldmax(alpha,f)

        if (efoldmax.lt.efold) cycle
        
        lnRhoRehMax = hni1_lnrhoreh_max(alpha,f,xend,Pstar)

        print *,'alpha=',alpha,'f/Mp=',f,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = hni1_x_star(alpha,f,xend,w,lnRhoReh,Pstar,bfoldstar)

           eps1 = hni1_epsilon_one(xstar,alpha,f)
           eps2 = hni1_epsilon_two(xstar,alpha,f)
           eps3 = hni1_epsilon_three(xstar,alpha,f)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('hni1_predic.dat',alpha,f,eps1,eps2,eps3,ns,r,Treh)

           call livewrite('hni1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/f,1._kp-alpha/))

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
  xEnd = hni1_x_endinf(alpha,f)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = hni1_x_rrad(alpha,f,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = hni1_epsilon_one(xstar,alpha,f)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  hni1_epsilon_one(xend,alpha,f)
     VendOverVstar = hni1_norm_potential(xend,alpha,f)/hni1_norm_potential(xstar,alpha,f)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = hni1_x_rreh(alpha,f,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = hni1_x_star(alpha,f,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo




end program hni1main
