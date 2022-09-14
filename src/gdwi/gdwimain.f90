!test the reheating derivation from slow-roll
program gdwimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use gdwisr, only : gdwi_epsilon_one, gdwi_epsilon_two, gdwi_epsilon_three
  use gdwireheat, only : gdwi_lnrhoreh_max, gdwi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use gdwisr, only : gdwi_norm_potential, gdwi_x_endinf
  use gdwireheat, only : gdwi_x_rreh, gdwi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,n
  integer :: npts = 20

  integer :: Nphi0
  ! real(kp) :: phi0min=2._kp*sqrt(2._kp)
  real(kp) :: phi0min=7.6_kp
  real(kp) :: phi0max=10._kp**3

  real(kp) :: phi0,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  
  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: p,V1
  real(kp) :: x, xmin,xmax

  
  Nphi0=26

  Pstar = powerAmpScalar


  n = 500
  xmin = 0._kp
  xmax = 1.5_kp

  p = 2
  phi0 = 1._kp
  
  do i=1,n
     x=xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)
     V1 = gdwi_norm_potential(x,p,phi0)

     call livewrite('gdwi_potential.dat',x,V1)

     eps1 = gdwi_epsilon_one(x,p,phi0)
     eps2 = gdwi_epsilon_two(x,p,phi0)
     eps3 = gdwi_epsilon_three(x,p,phi0)

     call livewrite('gdwi_slowroll.dat',x,eps1,eps2,eps3)

  end do

  
  call delete_file('gdwi_predic.dat')
  call delete_file('gdwi_nsr.dat')

  call aspicwrite_header('gdwi',labeps12,labnsr,labbfoldreh,(/'phi0','p   '/))
  
  !  w = 1._kp/3._kp
  w=0._kp

  p = 2._kp
  
  do j=0,Nphi0 

!     phi0=phi0min+(phi0max-phi0min)*(real(j,kp)/real(Nphi0,kp)) !arithmetic step
!     phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(Nphi0,kp)) !logarithmic step
!     phi0=exp(log(phi0min)*(log(phi0max)/log(phi0min))**(real(j,kp)/real(Nphi0,kp))) !superlogarithmic step
     phi0=exp(exp(log(log(phi0min))*(log(log(phi0max))/log(log(phi0min)))**(real(j,kp)/real(Nphi0,kp)))) !ultralogarithmic step



     lnRhoRehMin = lnRhoNuc
     xEnd = gdwi_x_endinf(p,phi0)
     lnRhoRehMax = gdwi_lnrhoreh_max(p,phi0,xend,Pstar)

     print *,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)



	xstar = gdwi_x_star(p,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)



        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar


        eps1 = gdwi_epsilon_one(xstar,p,phi0)
        eps2 = gdwi_epsilon_two(xstar,p,phi0)
        eps3 = gdwi_epsilon_three(xstar,p,phi0)


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)


        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('gdwi_predic.dat',phi0,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('gdwi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/phi0,p/))
        
     end do

  end do


   p = 3._kp
  
  do j=0,Nphi0 

!     phi0=phi0min+(phi0max-phi0min)*(real(j,kp)/real(Nphi0,kp)) !arithmetic step
!     phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(Nphi0,kp)) !logarithmic step
!     phi0=exp(log(phi0min)*(log(phi0max)/log(phi0min))**(real(j,kp)/real(Nphi0,kp))) !superlogarithmic step
     phi0=exp(exp(log(log(phi0min))*(log(log(phi0max))/log(log(phi0min)))**(real(j,kp)/real(Nphi0,kp)))) !ultralogarithmic step



     lnRhoRehMin = lnRhoNuc
     xEnd = gdwi_x_endinf(p,phi0)
     lnRhoRehMax = gdwi_lnrhoreh_max(p,phi0,xend,Pstar)

     print *,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)



	xstar = gdwi_x_star(p,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)



        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar


        eps1 = gdwi_epsilon_one(xstar,p,phi0)
        eps2 = gdwi_epsilon_two(xstar,p,phi0)
        eps3 = gdwi_epsilon_three(xstar,p,phi0)


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)


        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('gdwi_predic.dat',phi0,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('gdwi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/phi0,p/))
        
     end do

  end do


 p = 4._kp
  
  do j=0,Nphi0 

!     phi0=phi0min+(phi0max-phi0min)*(real(j,kp)/real(Nphi0,kp)) !arithmetic step
!     phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(Nphi0,kp)) !logarithmic step
!     phi0=exp(log(phi0min)*(log(phi0max)/log(phi0min))**(real(j,kp)/real(Nphi0,kp))) !superlogarithmic step
     phi0=exp(exp(log(log(phi0min))*(log(log(phi0max))/log(log(phi0min)))**(real(j,kp)/real(Nphi0,kp)))) !ultralogarithmic step



     lnRhoRehMin = lnRhoNuc
     xEnd = gdwi_x_endinf(p,phi0)
     lnRhoRehMax = gdwi_lnrhoreh_max(p,phi0,xend,Pstar)

     print *,'phi0= p=',phi0,p,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)



	xstar = gdwi_x_star(p,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)



        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar


        eps1 = gdwi_epsilon_one(xstar,p,phi0)
        eps2 = gdwi_epsilon_two(xstar,p,phi0)
        eps3 = gdwi_epsilon_three(xstar,p,phi0)


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)


        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('gdwi_predic.dat',phi0,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('gdwi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/phi0,p/))
        
     end do

  end do
  
  
  call aspicwrite_end()
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('gdwi_predic_summarized.dat') 
  nphi0=1000
  phi0min=7.6_kp
  phi0max=10._kp**2
  w=0._kp
  do j=1,nphi0
     phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(nphi0,kp))
     lnRhoReh = lnRhoNuc
     xEnd = gdwi_x_endinf(p,phi0)
     xstarA = gdwi_x_star(p,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)
     eps1A = gdwi_epsilon_one(xstarA,p,phi0)
     eps2A = gdwi_epsilon_two(xstarA,p,phi0)
     eps3A = gdwi_epsilon_three(xstarA,p,phi0)
     nsA = 1._kp - 2._kp*eps1A - eps2A
     rA = 16._kp*eps1A
     lnRhoReh = gdwi_lnrhoreh_max(p,phi0,xend,Pstar)
     xstarB = gdwi_x_star(p,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)
     eps1B = gdwi_epsilon_one(xstarB,p,phi0)
     eps2B = gdwi_epsilon_two(xstarB,p,phi0)
     eps3B = gdwi_epsilon_three(xstarB,p,phi0)
     nsB = 1._kp - 2._kp*eps1B - eps2B
     rB =16._kp*eps1B
     call livewrite('gdwi_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
  enddo


   write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  phi0 = 20

  xEnd = gdwi_x_endinf(p,phi0)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = gdwi_x_rrad(p,phi0,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = gdwi_epsilon_one(xstar,p,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  gdwi_epsilon_one(xend,p,phi0)
     VendOverVstar = gdwi_norm_potential(xend,p,phi0)/gdwi_norm_potential(xstar,p,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = gdwi_x_rreh(p,phi0,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = gdwi_x_star(p,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program gdwimain
