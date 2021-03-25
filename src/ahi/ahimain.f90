!test the reheating derivation from slow-roll
program ahimain
  use infprec, only : kp, pi
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use ahisr, only : ahi_epsilon_one, ahi_epsilon_two, ahi_epsilon_three
  use ahireheat, only : ahi_lnrhoreh_max, ahi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use ahisr, only : ahi_norm_potential, ahi_x_endinf
  use ahireheat, only : ahi_x_rreh, ahi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20,n

  real(kp) :: f,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(1:4) ::fvalues

  real(kp)  :: alpha,alphamin,alphamax,eps1A,eps2A,eps3A,nsA,rA
  real(kp)  :: eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  integer :: nalpha

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  
  real(kp) :: x,xmin,xmax,V1
  real(kp) :: fmin,fmax
  

  Pstar = powerAmpScalar

  
  call delete_file('ahi_potential.dat')
  call delete_file('ahi_slowroll.dat')

  n=250

  xmin = -2._kp
  xmax = 8.3_kp
  
  do i=1,n
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)        

     V1 = ahi_norm_potential(x,f=1._kp)


     call livewrite('ahi_potential.dat',x,V1)
     
     eps1 = ahi_epsilon_one(x,f=1._kp)
     eps2 = ahi_epsilon_two(x,f=1._kp)
     eps3 = ahi_epsilon_three(x,f=1._kp)
          
     call livewrite('ahi_slowroll.dat',x,eps1,eps2,eps3)
     
  enddo

  call delete_file('ahi_xend.dat')

  n=250
  fmin = 1d-6
  fmax = 1e3
  
  do i=1,n
     f=exp(log(fmin) + real(i-1,kp)*(log(fmax)-log(fmin))/real(n-1,kp))

     xend = ahi_x_endinf(f)

     call livewrite('ahi_xend.dat',f,xend)
     
  enddo
  
  
  call delete_file('ahi_predic.dat')
  call delete_file('ahi_nsr.dat')

  call aspicwrite_header('ahi',labeps12,labnsr,labbfoldreh,(/'f'/))

  fmin = 0.1
  fmax = 10.

  n = 20
  
  do j=1,n

     f=fmin + (fmax-fmin)*real(j-1,kp)/real(n-1,kp)

     w=0._kp

     lnRhoRehMin = lnRhoNuc

     xend = ahi_x_endinf(f)
     lnRhoRehMax = ahi_lnrhoreh_max(f,xend,Pstar)

     print *,'f=',f,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = ahi_x_star(f,xend,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = ahi_epsilon_one(xstar,f)
        eps2 = ahi_epsilon_two(xstar,f)
        eps3 = ahi_epsilon_three(xstar,f)


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(pi**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('ahi_predic.dat',f,eps1,eps2,eps3,r,ns,Treh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/f/))
        
     end do
  end do

  call aspicwrite_end()

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  f = 1.

  xend = ahi_x_endinf(f)
  
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = ahi_x_rrad(f,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = ahi_epsilon_one(xstar,f)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  ahi_epsilon_one(xend,f)
     VendOverVstar = ahi_norm_potential(xend,f)/ahi_norm_potential(xstar,f)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = ahi_x_rreh(f,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = ahi_x_star(f,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program ahimain
