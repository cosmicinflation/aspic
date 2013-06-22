!test the reheating derivation from slow-roll
program vhimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use vhisr, only : vhi_epsilon_one, vhi_epsilon_two, vhi_epsilon_three
  use vhisr, only : vhi_xinimax,vhi_xendmin,vhi_xendmax
  use vhireheat, only : vhi_lnrhoreh_max, vhi_x_star
  use infinout, only : delete_file, livewrite, has_not_shifted
  use srreheat, only : log_energy_reheat_ingev

  use vhisr, only : vhi_norm_potential
  use vhireheat, only : vhi_x_rreh, vhi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 20

  integer :: Nmu
  real(kp) :: mumin
  real(kp) :: mumax

  integer :: NxEnd
  real(kp) :: xEndmin
  real(kp) :: xEndmax

  real(kp) :: p,mu,xEnd,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End

  Pstar = powerAmpScalar

  call delete_file('vhi_predic.dat')
  call delete_file('vhi_nsr.dat')


  !  w = 1._kp/3._kp
  w=0._kp




!!!!!!!!!!!!!!!!!!!!!
!!!!    p=0.5    !!!!
!!!!!!!!!!!!!!!!!!!!!

  p=0.5_kp
  Nmu=8
  mumin=0.001_kp
  mumax=1000000._kp
  NxEnd=30

  do j=0,Nmu
     mu=mumin*(mumax/mumin)**(real(j,kp)/Nmu)
     xEndmin=vhi_xendmin(p,mu)*(1._kp+epsilon(1._kp))
     xEndmax=vhi_xinimax(p,mu)*0.99

     do k=0,NxEnd
        xEnd=xEndmin*(xEndmax/xEndmin)**(real(k,kp)/NxEnd) !logarithmic step
        ! xEnd=xEndmin+(xEndmax-xEndmin)*(real(k,kp)/NxEnd)*p,mu*(1.+epsilon(1._kp)) !arithmetic step


        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = vhi_lnrhoreh_max(p,mu,xEnd,Pstar)

        print *,'p',p,'mu=',mu,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = vhi_x_star(p,mu,xEnd,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',vhi_epsilon_one(xstar,p,mu)



           eps1 = vhi_epsilon_one(xstar,p,mu)
           eps2 = vhi_epsilon_two(xstar,p,mu)
           eps3 = vhi_epsilon_three(xstar,p,mu)


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)


           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

            if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
              cycle
           endif

           call livewrite('vhi_predic.dat',p,mu,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('vhi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do

!!!!!!!!!!!!!!!!!!!!!
!!!!     p=1     !!!!
!!!!!!!!!!!!!!!!!!!!!


  p=1._kp
  Nmu=8
  mumin=0.001_kp
  mumax=1000._kp
  NxEnd=15

  do j=0,Nmu
     mu=mumin*(mumax/mumin)**(real(j,kp)/Nmu)
     xEndmin=vhi_xendmin(p,mu)*(1._kp+epsilon(1._kp))
     xEndmax=vhi_xinimax(p,mu)*0.99

     do k=0,NxEnd
        xEnd=xEndmin*(xEndmax/xEndmin)**(real(k,kp)/NxEnd) !logarithmic step
        ! xEnd=xEndmin+(xEndmax-xEndmin)*(real(k,kp)/NxEnd)*p,mu*(1.+epsilon(1._kp)) !arithmetic step


        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = vhi_lnrhoreh_max(p,mu,xEnd,Pstar)

        print *,'p',p,'mu=',mu,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = vhi_x_star(p,mu,xEnd,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',vhi_epsilon_one(xstar,p,mu)



           eps1 = vhi_epsilon_one(xstar,p,mu)
           eps2 = vhi_epsilon_two(xstar,p,mu)
           eps3 = vhi_epsilon_three(xstar,p,mu)


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)


           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

            if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
              cycle
           endif

           call livewrite('vhi_predic.dat',p,mu,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('vhi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do


!!!!!!!!!!!!!!!!!!!!!
!!!!    p=1.5    !!!!
!!!!!!!!!!!!!!!!!!!!!

  p=1.5_kp
  Nmu=15
  mumin=0.8_kp
  mumax=1000._kp
  NxEnd=40

  do j=0,Nmu
     mu=mumin*(mumax/mumin)**(real(j,kp)/Nmu)
     xEndmin=vhi_xendmin(p,mu)*(1._kp+epsilon(1._kp))
     xEndmax=vhi_xinimax(p,mu)*0.99
     xEndmax=min(vhi_xendmax(50._kp,p,mu),10._kp) !To remain vacum dominated

     do k=0,NxEnd
        xEnd=xEndmin*(xEndmax/xEndmin)**(real(k,kp)/NxEnd) !logarithmic step
        ! xEnd=xEndmin+(xEndmax-xEndmin)*(real(k,kp)/NxEnd)*p,mu*(1.+epsilon(1._kp)) !arithmetic step


        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = vhi_lnrhoreh_max(p,mu,xEnd,Pstar)

        print *,'p,mu=',p,mu,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = vhi_x_star(p,mu,xEnd,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',vhi_epsilon_one(xstar,p,mu)



           eps1 = vhi_epsilon_one(xstar,p,mu)
           eps2 = vhi_epsilon_two(xstar,p,mu)
           eps3 = vhi_epsilon_three(xstar,p,mu)

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)


           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
              cycle
           endif

           call livewrite('vhi_predic.dat',p,mu,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('vhi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do


!!!!!!!!!!!!!!!!!!!!!
!!!!     p=2     !!!!
!!!!!!!!!!!!!!!!!!!!!

  p=2._kp
  Nmu=25
  mumin=0.8_kp
  mumax=1000._kp
  NxEnd=50

  do j=0,Nmu
     mu=mumin*(mumax/mumin)**(real(j,kp)/Nmu)
     xEndmin=vhi_xendmin(p,mu)*(1._kp+epsilon(1._kp))
     xEndmax=vhi_xinimax(p,mu)*0.99
     xEndmax=min(vhi_xendmax(50._kp,p,mu),10._kp) !To remain vacum dominated

     do k=0,NxEnd
        xEnd=xEndmin*(xEndmax/xEndmin)**(real(k,kp)/NxEnd) !logarithmic step
        ! xEnd=xEndmin+(xEndmax-xEndmin)*(real(k,kp)/NxEnd)*p,mu*(1.+epsilon(1._kp)) !arithmetic step


        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = vhi_lnrhoreh_max(p,mu,xEnd,Pstar)

        print *,'p,mu=',p,mu,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = vhi_x_star(p,mu,xEnd,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',vhi_epsilon_one(xstar,p,mu)



           eps1 = vhi_epsilon_one(xstar,p,mu)
           eps2 = vhi_epsilon_two(xstar,p,mu)
           eps3 = vhi_epsilon_three(xstar,p,mu)


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)


           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
              cycle
           endif

           call livewrite('vhi_predic.dat',p,mu,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('vhi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do


!!!!!!!!!!!!!!!!!!!!!
!!!!     p=3     !!!!
!!!!!!!!!!!!!!!!!!!!!

  p=3._kp
  Nmu=20
  mumin=0.8_kp
  mumax=1000._kp
  NxEnd=200

  do j=0,Nmu
     mu=mumin*(mumax/mumin)**(real(j,kp)/Nmu)
     xEndmin=vhi_xendmin(p,mu)*(1._kp+epsilon(1._kp))
     xEndmax=vhi_xinimax(p,mu)*0.99
     xEndmax=min(vhi_xendmax(50._kp,p,mu),10._kp) !To remain vacum dominated

     do k=0,NxEnd
        xEnd=xEndmin*(xEndmax/xEndmin)**(real(k,kp)/NxEnd) !logarithmic step
        ! xEnd=xEndmin+(xEndmax-xEndmin)*(real(k,kp)/NxEnd)*p,mu*(1.+epsilon(1._kp)) !arithmetic step


        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = vhi_lnrhoreh_max(p,mu,xEnd,Pstar)

        print *,'p,mu=',p,mu,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = vhi_x_star(p,mu,xEnd,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',vhi_epsilon_one(xstar,p,mu)



           eps1 = vhi_epsilon_one(xstar,p,mu)
           eps2 = vhi_epsilon_two(xstar,p,mu)
           eps3 = vhi_epsilon_three(xstar,p,mu)


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)


           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )


           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
              cycle
           endif

           call livewrite('vhi_predic.dat',p,mu,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('vhi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do





  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  p=3.
  mu =100.
  xend =10.
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = vhi_x_rrad(p,mu,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = vhi_epsilon_one(xstar,p,mu)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  vhi_epsilon_one(xend,p,mu)
     VendOverVstar = vhi_norm_potential(xend,p,mu)/vhi_norm_potential(xstar,p,mu)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = vhi_x_rreh(p,mu,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = vhi_x_star(p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program vhimain
