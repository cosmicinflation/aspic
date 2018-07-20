!test the reheating derivation from slow-roll
program rpi3main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rpi3sr, only : rpi3_epsilon_one, rpi3_epsilon_two, rpi3_epsilon_three
  use rpi3reheat, only : rpi3_lnrhoreh_max, rpi3_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use rpi3sr, only : rpi3_norm_potential, rpi3_x_endinf
  use rpi3reheat, only : rpi3_x_rreh, rpi3_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 15

  integer :: Np=20
  real(kp) :: pmin=0.8_kp
  real(kp) :: pmax=1._kp

  real(kp) :: p,w,bfoldstar
  real(kp) :: lnRhoReh,ystar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar

  call delete_file('rpi3_predic.dat')
  call delete_file('rpi3_nsr.dat')

  call aspicwrite_header('rpi3',labeps12,labnsr,labbfoldreh,(/'1mp'/))
  
  !  w = 1._kp/3._kp
  w=0._kp

  do j=0,Np 
     p=pmin+(pmax-pmin)*(real(j,kp)/Np)
!     w=(1._kp-p)/(3._kp-1._kp)

     print *,'p= ',p

     lnRhoRehMin = lnRhoNuc
     xEnd = rpi3_x_endinf(p)
     lnRhoRehMax = rpi3_lnrhoreh_max(p,xend,Pstar)

     print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)



        ystar = rpi3_x_star(p,xend,w,lnRhoReh,Pstar,bfoldstar)


        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'ystar=',ystar


        eps1 = rpi3_epsilon_one(ystar,p)
        eps2 = rpi3_epsilon_two(ystar,p)
        eps3 = rpi3_epsilon_three(ystar,p)


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)


        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/1._kp-p/))
        
        call livewrite('rpi3_predic.dat',p,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('rpi3_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do

  call aspicwrite_end()
  
  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  p = 0.9
  xEnd = rpi3_x_endinf(p)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     ystar = rpi3_x_rrad(p,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'ystar', ystar

     eps1 = rpi3_epsilon_one(ystar,p)

     !consistency test
     !get lnR from lnRrad and check that it gives the same ystar
     eps1end =  rpi3_epsilon_one(xend,p)
     VendOverVstar = rpi3_norm_potential(xend,p)/rpi3_norm_potential(ystar,p)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     ystar = rpi3_x_rreh(p,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'ystar', ystar

     !second consistency check
     !get rhoreh for chosen w and check that ystar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     ystar = rpi3_x_star(p,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'ystar',ystar

  enddo

end program rpi3main
