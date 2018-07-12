!test the reheating derivation from slow-roll
program hbimain
  use infprec, only : kp, transfert
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use hbireheat, only : hbi_x_star, hbi_lnrhoreh_max
  use hbireheat, only : hbi_xnphi0_fromepsilon, hbi_lnrhoreh_fromepsilon
  use hbisr, only : hbi_epsilon_one, hbi_epsilon_two,hbi_epsilon_three
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use hbisr, only : hbi_norm_potential, hbi_x_endinf
  use hbireheat, only : hbi_x_rreh, hbi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat
  
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none

  type(transfert) :: hbiData
  real(kp) :: Pstar,calF

  integer :: i,j
  integer :: npts = 20,nphi0=30.

  real(kp) :: phi0,n,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r,Treh
  real(kp) :: logErehGeV

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(3) :: vecbuffer

  real(kp) ::phi0min,phi0max

  logical, parameter :: display = .true.
  logical, parameter :: inversion = .true.

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  w = 0._kp

  Pstar = powerAmpScalar

  call delete_file('hbi_predic.dat')
  call delete_file('hbi_nsr.dat')

  call aspicwrite_header('hbi',labeps12,labnsr,labbfoldreh,(/'phi0','n   '/))
  
!!!!!!!!!!!!!! 
!!!! n=1  !!!!
!!!!!!!!!!!!!!

  n = 1._kp

  phi0min=n/sqrt(2.)*10._kp
  phi0max=10._kp**(3.)

  do j=0,nphi0

     ! Ultralogarithmic step
     phi0=phi0min/exp(1._kp)*exp(exp(real(j,kp)/real(nphi0,kp)*log(1._kp+log(phi0max/phi0min))))


     !xstar stands for phistar/phi0

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = hbi_lnrhoreh_max(n,phi0,Pstar)

     print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = hbi_x_star(n,phi0,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = hbi_epsilon_one(xstar,n,phi0)
        eps2 = hbi_epsilon_two(xstar,n,phi0)
        eps3 = hbi_epsilon_three(xstar,n,phi0)

        if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfoldstar),eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('hbi_predic.dat',n,phi0,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('hbi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/phi0,n/))
        
     end do

  end do

!!!!!!!!!!!!!! 
!!!! n=2  !!!!
!!!!!!!!!!!!!!

  n = 2._kp

  phi0min=n/sqrt(2.)*10._kp
  phi0max=10._kp**(3.)

  do j=0,nphi0

     ! Ultralogarithmic step
     phi0=phi0min/exp(1._kp)*exp(exp(real(j,kp)/real(nphi0,kp)*log(1._kp+log(phi0max/phi0min))))


     !xstar stands for phistar/phi0

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = hbi_lnrhoreh_max(n,phi0,Pstar)

     print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = hbi_x_star(n,phi0,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = hbi_epsilon_one(xstar,n,phi0)
        eps2 = hbi_epsilon_two(xstar,n,phi0)
        eps3 = hbi_epsilon_three(xstar,n,phi0)

        if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfoldstar),eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('hbi_predic.dat',n,phi0,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('hbi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/phi0,n/))

     end do

  end do

!!!!!!!!!!!!!! 
!!!! n=3  !!!!
!!!!!!!!!!!!!!

  n = 3._kp

  phi0min=n/sqrt(2.)*10._kp
  phi0max=10._kp**(3.)

  do j=0,nphi0

     ! Ultralogarithmic step
     phi0=phi0min/exp(1._kp)*exp(exp(real(j,kp)/real(nphi0,kp)*log(1._kp+log(phi0max/phi0min))))


     !xstar stands for phistar/phi0

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = hbi_lnrhoreh_max(n,phi0,Pstar)

     print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = hbi_x_star(n,phi0,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = hbi_epsilon_one(xstar,n,phi0)
        eps2 = hbi_epsilon_two(xstar,n,phi0)
        eps3 = hbi_epsilon_three(xstar,n,phi0)

        if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfoldstar),eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('hbi_predic.dat',n,phi0,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('hbi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/phi0,n/))
        
     end do

  end do




  call aspicwrite_end()





  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  n=2.5
  phi0 = 10._kp
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = hbi_x_rrad(n,phi0,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = hbi_epsilon_one(xstar,n,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = hbi_x_endinf(n,phi0)
     eps1end =  hbi_epsilon_one(xend,n,phi0)
     VendOverVstar = hbi_norm_potential(xend,n,phi0)/hbi_norm_potential(xstar,n,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = hbi_x_rreh(n,phi0,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = hbi_x_star(n,phi0,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program hbimain
