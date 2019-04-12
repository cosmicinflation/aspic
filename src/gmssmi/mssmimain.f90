!test the reheating derivation from slow-roll
program mssmimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use mssmisr, only : mssmi_epsilon_one, mssmi_epsilon_two, mssmi_epsilon_three, mssmi_x_epsonemin, mssmi_x_endinf
  use mssmireheat, only : mssmi_lnrhoreh_max, mssmi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use mssmisr, only : mssmi_norm_potential, mssmi_x_endinf
  use mssmireheat, only : mssmi_x_rreh, mssmi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts,nphi0

  real(kp) :: phi0,w,bfoldstar,phi0min,phi0max
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::x,xmin,xmax,xendNUM,xendANAL


  real(kp) :: eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        Tests the approximated formula for xend        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  phi0min=10._kp**(-6.)
  phi0max=10._kp**(-0.)
  nphi0=100
  call delete_file('mssmi_xend.dat')

  do j=0,nphi0
     phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(nphi0,kp))
     xendNUM=mssmi_x_endinf(phi0)
     xendANAL=1._kp-2._kp**(-0.75_kp)*sqrt(phi0/15._kp)
     call livewrite('mssmi_xend.dat',phi0,xendNUM,xendANAL)
  end do
  print*,'mssmi_xend.dat written.'



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts = 20
  nphi0=100

  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('mssmi_predic.dat')
  call delete_file('mssmi_nsr.dat')

  call aspicwrite_header('mssmi',labeps12,labnsr,labbfoldreh,(/'phi0'/))
  
  !Prior on phi0
  phi0min=10._kp**(-2.)
  phi0max=10._kp**(3.)

  call livewrite('mssmi_predic.dat',10._kp**(-4.), &
       10._kp**(-8)/(60._kp*50._kp**2),4._kp/50._kp, &
       1._kp/50._kp,16._kp*10._kp**(-8)/(60._kp*50._kp**2), &
       1._kp -4._kp/50._kp,1._kp) !To prime the color bar

  do j=0,nphi0
     phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(nphi0,kp))


     lnRhoRehMin = lnRhoNuc
     xEnd = mssmi_x_endinf(phi0)
     lnRhoRehMax = mssmi_lnrhoreh_max(phi0,xend,Pstar)

     print *,'phi0=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = mssmi_x_star(phi0,xend,w,lnRhoReh,Pstar,bfoldstar)


        eps1 = mssmi_epsilon_one(xstar,phi0)
        eps2 = mssmi_epsilon_two(xstar,phi0)
        eps3 = mssmi_epsilon_three(xstar,phi0)


        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('mssmi_predic.dat',phi0,eps1,eps2,eps3,r,ns,Treh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/phi0/))

     end do


  end do


  call aspicwrite_end()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('mssmi_predic_summarized.dat') 
  nphi0=1000
  phi0min=10._kp**(-2.)
  phi0max=10._kp**(3.)
  w=0._kp
  do j=0,nphi0
     phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(nphi0,kp))
     xEnd = mssmi_x_endinf(phi0)
     lnRhoReh = lnRhoNuc
     xstarA = mssmi_x_star(phi0,xend,w,lnRhoReh,Pstar,bfoldstar)
     eps1A = mssmi_epsilon_one(xstarA,phi0)
     eps2A = mssmi_epsilon_two(xstarA,phi0)
     eps3A = mssmi_epsilon_three(xstarA,phi0)
     nsA = 1._kp - 2._kp*eps1A - eps2A
     rA = 16._kp*eps1A
     lnRhoReh = mssmi_lnrhoreh_max(phi0,xend,Pstar)
     xstarB = mssmi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)
     eps1B = mssmi_epsilon_one(xstarB,phi0)
     eps2B = mssmi_epsilon_two(xstarB,phi0)
     eps3B = mssmi_epsilon_three(xstarB,phi0)
     nsB = 1._kp - 2._kp*eps1B - eps2B
     rB =16._kp*eps1B
     if ((rA .gt. 0._kp) .and. (nsA .gt.0._kp) .and. &
          (rB .gt. 0._kp) .and. (nsB .gt. 0._kp) &
          .and. (eps1A .gt. 10.**(-10.)) .and. (eps1B .gt. 10.**(-10.)) &
          .and. (eps2A .gt. 10.**(-10.)) .and. (eps2B .gt. 10.**(-10.))) then !to remove NaN and error points
        call livewrite('mssmi_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
     end if
  enddo


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  phi0 = 10
  xEnd = mssmi_x_endinf(phi0)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = mssmi_x_rrad(phi0,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = mssmi_epsilon_one(xstar,phi0)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  mssmi_epsilon_one(xend,phi0)
     VendOverVstar = mssmi_norm_potential(xend,phi0)/mssmi_norm_potential(xstar,phi0)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = mssmi_x_rreh(phi0,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = mssmi_x_star(phi0,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program mssmimain
