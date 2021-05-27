!test the reheating derivation from slow-roll
program simain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use sisr, only : si_epsilon_one, si_epsilon_two, si_epsilon_three
  use sireheat, only : si_lnrhoreh_max, si_x_star
  use infinout, only : delete_file, livewrite
  use infinout, only : labeps12, labnsr, labbfoldreh
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use srreheat, only : log_energy_reheat_ingev

  use sisr, only : si_norm_potential, si_x_endinf
  use sireheat, only : si_x_rreh, si_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,n
  integer :: npts = 20

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend
  real(kp) :: lnOmega4End

  real(kp) :: xmin, xmax, V1,x

  n = 100

  call delete_file('si_potential.dat')
  call delete_file('si_slowroll.dat')

  
  n=250

  xmin = 0.0_kp
  xmax = 10._kp
  
  do i=1,n
     
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)

     V1 = si_norm_potential(x)

     call livewrite('si_potential.dat',x,V1)
     
     eps1 = si_epsilon_one(x)
     eps2 = si_epsilon_two(x)
     eps3 = si_epsilon_three(x)

          
     call livewrite('si_slowroll.dat',x,eps1,eps2,eps3)
     
  enddo
  
  Pstar = powerAmpScalar

  call delete_file('si_predic.dat')
  call delete_file('si_nsr.dat')


  !  w = 1._kp/3._kp
  w=0._kp

  call aspicwrite_header('si',labeps12,labnsr,labbfoldreh)
  
  lnRhoRehMin = lnRhoNuc
  xEnd = si_x_endinf()
  lnRhoRehMax = si_lnrhoreh_max(xend,Pstar)

  print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)



     xstar = si_x_star(xend,w,lnRhoReh,Pstar,bfoldstar)



     print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar


     eps1 = si_epsilon_one(xstar)
     eps2 = si_epsilon_two(xstar)
     eps3 = si_epsilon_three(xstar)


     logErehGeV = log_energy_reheat_ingev(lnRhoReh)


     Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


     ns = 1._kp - 2._kp*eps1 - eps2
     r =16._kp*eps1

     call livewrite('si_predic.dat',eps1,eps2,eps3,r,ns,Treh)

     call livewrite('si_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/))
     
  end do

  call aspicwrite_end()

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('si_predic_summarized.dat')
  lnRhoReh = lnRhoRehMin
  xEnd = si_x_endinf()
  xstar = si_x_star(xend,w,lnRhoReh,Pstar,bfoldstar)
  eps1 = si_epsilon_one(xstar)
  eps2 = si_epsilon_two(xstar)
  eps3 = si_epsilon_three(xstar)
  ns = 1._kp - 2._kp*eps1 - eps2
  r =16._kp*eps1
  call livewrite('si_predic_summarized.dat',eps1,eps2,eps3,r,ns)
  lnRhoReh = lnRhoRehMax
  xstar = si_x_star(xend,w,lnRhoReh,Pstar,bfoldstar)
  eps1 = si_epsilon_one(xstar)
  eps2 = si_epsilon_two(xstar)
  eps3 = si_epsilon_three(xstar)
  ns = 1._kp - 2._kp*eps1 - eps2
  r =16._kp*eps1
  call livewrite('si_predic_summarized.dat',eps1,eps2,eps3,r,ns)

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  xEnd = si_x_endinf()
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = si_x_rrad(xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = si_epsilon_one(xstar)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  si_epsilon_one(xend)
     VendOverVstar = si_norm_potential(xend)/si_norm_potential(xstar)
     lnOmega4End =  2._kp*sqrt(2._kp/3._kp)*xend

!Jordan frame (we input the conformal factor)     
     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar,lnOmega4End)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = si_x_rreh(xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
!Jordan frame (we input the conformal factor)
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar,lnOmega4End)

     xstar = si_x_star(xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program simain
