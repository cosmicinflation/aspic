!test the reheating derivation from slow-roll
program hbimain
  use infprec, only : kp, transfert
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use hbireheat, only : hbi_x_star, hbi_lnrhoreh_max
  use hbireheat, only : hbi_xnmu_fromepsilon, hbi_lnrhoreh_fromepsilon
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
  integer :: npts = 20,nmu=30.

  real(kp) :: mu,n,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r,Treh
  real(kp) :: logErehGeV

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(3) :: vecbuffer

  real(kp) ::mumin,mumax

  logical, parameter :: display = .true.
  logical, parameter :: inversion = .true.

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: epsOneEnd,potEnd ! Used for bebugging

  integer :: npt
  real(kp) :: xmin, xmax, x,V1
  
  w = 0._kp

  Pstar = powerAmpScalar

  call delete_file('hbi_potential.dat')
  call delete_file('hbi_slowroll.dat')



  npt=250

  xmin = 0._kp
  xmax = 5._kp

  do i=1,npt
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(npt-1,kp)        

     V1 = hbi_norm_potential(x,n=1.1_kp,mu=1._kp)


     call livewrite('hbi_potential.dat',x,V1)

     eps1 = hbi_epsilon_one(x,n=1.1_kp,mu=1._kp)
     eps2 = hbi_epsilon_two(x,n=1.1_kp,mu=1._kp)
     eps3 = hbi_epsilon_three(x,n=1.1_kp,mu=1._kp)

     call livewrite('hbi_slowroll.dat',x,eps1,eps2,eps3)


  enddo



  

  call delete_file('hbi_predic.dat')
  call delete_file('hbi_nsr.dat')

  call aspicwrite_header('hbi',labeps12,labnsr,labbfoldreh,(/'mu  ','n   '/))
  
!!!!!!!!!!!!!!!! 
!!!! n=0.5  !!!!
!!!!!!!!!!!!!!!!

  n = 0.5_kp

  mumin=n/sqrt(2.)*10._kp
  mumax=100._kp

  do j=0,nmu

     ! Ultralogarithmic step
     mu=mumin/exp(1._kp)*exp(exp(real(j,kp)/real(nmu,kp)*log(1._kp+log(mumax/mumin))))


     !xstar stands for phistar/mu

     lnRhoRehMin = lnRhoNuc
     xEnd = hbi_x_endinf(n,mu)       
     lnRhoRehMax = hbi_lnrhoreh_max(n,mu,xend,Pstar)

     print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

     do i=npts,1,-1

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = hbi_x_star(n,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = hbi_epsilon_one(xstar,n,mu)
        eps2 = hbi_epsilon_two(xstar,n,mu)
        eps3 = hbi_epsilon_three(xstar,n,mu)

        if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfoldstar),eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('hbi_predic.dat',n,mu,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('hbi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/mu,n/))
        
     end do

  end do

!!!!!!!!!!!!!! 
!!!! n=1  !!!!
!!!!!!!!!!!!!!

  n = 1._kp

  mumin=n/sqrt(2.)*10._kp
  mumax=10._kp**(2.)

  do j=0,nmu

     ! Ultralogarithmic step
     mu=mumin/exp(1._kp)*exp(exp(real(j,kp)/real(nmu,kp)*log(1._kp+log(mumax/mumin))))


     !xstar stands for phistar/mu

     lnRhoRehMin = lnRhoNuc
     xEnd = hbi_x_endinf(n,mu)       
     lnRhoRehMax = hbi_lnrhoreh_max(n,mu,xend,Pstar)

     print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

     do i=npts,1,-1

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = hbi_x_star(n,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = hbi_epsilon_one(xstar,n,mu)
        eps2 = hbi_epsilon_two(xstar,n,mu)
        eps3 = hbi_epsilon_three(xstar,n,mu)

        if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfoldstar),eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('hbi_predic.dat',n,mu,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('hbi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/mu,n/))

     end do

  end do

!!!!!!!!!!!!!!!! 
!!!! n=1.5  !!!!
!!!!!!!!!!!!!!!!

  n = 1.5_kp

  mumin=n/sqrt(2.)*10._kp
  mumax=10._kp**(2.)

  do j=0,nmu

     ! Ultralogarithmic step
     mu=mumin/exp(1._kp)*exp(exp(real(j,kp)/real(nmu,kp)*log(1._kp+log(mumax/mumin))))


     !xstar stands for phistar/mu

     lnRhoRehMin = lnRhoNuc
     xEnd = hbi_x_endinf(n,mu)       
     lnRhoRehMax = hbi_lnrhoreh_max(n,mu,xend,Pstar)

     print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

     do i=npts,1,-1

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = hbi_x_star(n,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
        eps1 = hbi_epsilon_one(xstar,n,mu)
        eps2 = hbi_epsilon_two(xstar,n,mu)
        eps3 = hbi_epsilon_three(xstar,n,mu)

        if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfoldstar),eps1


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('hbi_predic.dat',n,mu,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('hbi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/mu,n/))
        
     end do

  end do




  call aspicwrite_end()





  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  n=2.5
  mu = 10._kp
  xEnd = hbi_x_endinf(n,mu)
  
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = hbi_x_rrad(n,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = hbi_epsilon_one(xstar,n,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  hbi_epsilon_one(xend,n,mu)
     VendOverVstar = hbi_norm_potential(xend,n,mu)/hbi_norm_potential(xstar,n,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = hbi_x_rreh(n,mu,xend,lnR,bfoldstar)
     !print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = hbi_x_star(n,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program hbimain
