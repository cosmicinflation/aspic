!test the reheating derivation from slow-roll
program hni2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar

  use hni2sr, only : hni2_norm_potential, hni2_numacc_efoldmax
  use hni2sr, only : hni2_epsilon_one, hni2_epsilon_two,hni2_epsilon_three
  use hni2sr, only : hni2_alphamax, hni2_fmin, hni2_numacc_xendmin, hni2_numacc_xendmax
  
  use hni2reheat, only : hni2_x_rreh, hni2_x_rrad
  use hni2reheat, only : hni2_lnrhoreh_max, hni2_x_star

  use hnicommon, only : hni_x_potmin
  
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat
  use srreheat, only : log_energy_reheat_ingev

  use infinout, only : delete_file, livewrite
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh

  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,l
  integer :: npts,nalpha,nf,ne

  real(kp) :: alpha,f,w,bfoldstar,alphamin,alphamax
  real(kp) :: fmin,fmax,fcutmin,fcutmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(:), allocatable ::fvalues

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: xendmin, xendmax, efold, efoldmax
  
  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  efold = 120._kp
  
  npts = 15
  nalpha= 3
  nf = 3
  
  ne = 100
  
  fcutmin = 0.2_kp
  fcutmax = 12._kp
  
  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('hni2_predic.dat')
  call delete_file('hni2_nsr.dat')

  call aspicwrite_header('hni2',labeps12,labnsr,labbfoldreh,(/'xend ','f    ','alpha'/))



  alphamin = 0.1
  alphamax = 0.97


  fmax= fcutmax
  fmin = fcutmin

  
  do k=1,nalpha
     alpha= alphamin + (alphamax - alphamin)*real(k-1,kp)/real(nalpha-1,kp)

     do j=1,nf

        fmin = max(2.0*hni2_fmin(alpha),fcutmin)
        
        f=fmin*(fmax/fmin)**(real(j-1,kp)/real(nf-1,kp))

        
        
        lnRhoRehMin = lnRhoNuc
        xendmax = 0.999*hni2_numacc_xendmax(alpha,f)

        efoldmax = hni2_numacc_efoldmax(alpha,f)

        if (efoldmax.lt.efold) cycle

        xendmin = hni2_numacc_xendmin(efold,alpha,f)
        
        
        do l=1,ne

           xend = xendmin + (xendmax-xendmin)*real(l-1,kp)/real(Ne-1,kp)
           
           lnRhoRehMax = hni2_lnrhoreh_max(alpha,f,xend,Pstar)

           print *,'alpha=',alpha,'f/Mp=',f,'xend',xend,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax


           do i=1,npts

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

              xstar = hni2_x_star(alpha,f,xend,w,lnRhoReh,Pstar,bfoldstar)

              eps1 = hni2_epsilon_one(xstar,alpha,f)
              eps2 = hni2_epsilon_two(xstar,alpha,f)
              eps3 = hni2_epsilon_three(xstar,alpha,f)

              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1


              logErehGeV = log_energy_reheat_ingev(lnRhoReh)
              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              call livewrite('hni2_predic.dat',alpha,f,eps1,eps2,eps3,ns,r,Treh)

              call livewrite('hni2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

              call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend,f,alpha/))

           end do

        enddo

     end do

  end do

  
  call aspicwrite_end()

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 0.05
  f = 10.
  xend = 2._kp
  
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = hni2_x_rrad(alpha,f,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = hni2_epsilon_one(xstar,alpha,f)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  hni2_epsilon_one(xend,alpha,f)
     VendOverVstar = hni2_norm_potential(xend,alpha,f)/hni2_norm_potential(xstar,alpha,f)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = hni2_x_rreh(alpha,f,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = hni2_x_star(alpha,f,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo




end program hni2main
