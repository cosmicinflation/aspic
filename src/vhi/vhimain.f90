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

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts = 30

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

  integer, parameter :: nvec = 2
  real(kp), dimension(nvec) :: muvec
  
  Pstar = powerAmpScalar

  
  call delete_file('vhi_predic.dat')
  call delete_file('vhi_nsr.dat')


  !  w = 1._kp/3._kp
  w=0._kp


  call aspicwrite_header('vhihalf',labeps12,labnsr,labbfoldreh,(/'xend','mu  ','p   '/))

!!!!!!!!!!!!!!!!!!!!!
!!!!    p=0.5    !!!!
!!!!!!!!!!!!!!!!!!!!!

  p=0.5_kp
  NxEnd=100

!  muvec = (/0.01, 10.0, 100.0/)

    muvec = (/0.01, 20.0/)
  
  do j=1,nvec
     mu=muvec(j)
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

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend,mu,p/))
           
            if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
              cycle
           endif

           call livewrite('vhi_predic.dat',p,mu,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('vhi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do

  call aspicwrite_end()
  
!!!!!!!!!!!!!!!!!!!!!
!!!!     p=1     !!!!
!!!!!!!!!!!!!!!!!!!!!

  call aspicwrite_header('vhione',labeps12,labnsr,labbfoldreh,(/'xend','mu  ','p   '/))
  
  p=1._kp
  NxEnd=100

  muvec = (/0.1, 2.0/)

  do j=1,nvec
     mu = muvec(j)
     xEndmin=vhi_xendmin(p,mu)*(1._kp+epsilon(1._kp))
     xEndmax=vhi_xinimax(p,mu)*0.1

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

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend,mu,p/))
           
            if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
              cycle
           endif

           call livewrite('vhi_predic.dat',p,mu,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('vhi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do

  call aspicwrite_end()

  
!!!!!!!!!!!!!!!!!!!!!
!!!!    p=1.5    !!!!
!!!!!!!!!!!!!!!!!!!!!

  call aspicwrite_header('vhithreehalf',labeps12,labnsr,labbfoldreh,(/'xend','mu  ','p   '/))
  
  p=1.5_kp
  NxEnd=100

!  muvec = (/1.0, 5.0, 7.0/)
  muvec = (/1.0, 5.0/)
  do j=1,nvec
     mu=muvec(j)

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

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend,mu,p/))
           
           if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
              cycle
           endif

           call livewrite('vhi_predic.dat',p,mu,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('vhi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do

  call aspicwrite_end()

!!!!!!!!!!!!!!!!!!!!!
!!!!     p=2     !!!!
!!!!!!!!!!!!!!!!!!!!!

  call aspicwrite_header('vhitwo',labeps12,labnsr,labbfoldreh,(/'xend','mu  ','p   '/))
  
  p=2._kp
  NxEnd=100

!  muvec = (/0.8, 8.0, 80.0/)

  muvec = (/0.8, 8.0/)

  
  do j=1,nvec
     mu=muvec(j)
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

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend,mu,p/))
           
           if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
              cycle
           endif

           if ((abs(eps2).gt.0.2).or.(eps1.lt.1e-6)) cycle

           call livewrite('vhi_predic.dat',p,mu,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('vhi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do

  call aspicwrite_end()

!!!!!!!!!!!!!!!!!!!!!
!!!!     p=3     !!!!
!!!!!!!!!!!!!!!!!!!!!

  call aspicwrite_header('vhithree',labeps12,labnsr,labbfoldreh,(/'xend','mu  ','p   '/))
  
  p=3._kp
  NxEnd=300

  muvec = (/5.0, 10.0/)
  
  do j=1,nvec
     mu=muvec(j)
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

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend,mu,p/))
           
           if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
              cycle
           endif

           if ((abs(eps2).gt.0.2).or.(eps1.lt.1e-6)) cycle

           call livewrite('vhi_predic.dat',p,mu,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('vhi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

     end do

  end do


  call aspicwrite_end()


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
