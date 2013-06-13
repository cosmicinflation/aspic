!test the reheating derivation from slow-roll
program cndimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use cndisr, only : cndi_epsilon_one, cndi_epsilon_two, cndi_epsilon_three, cndi_xendmax
  use cndireheat, only : cndi_lnrhoreh_max, cndi_x_star
  use infinout, only : delete_file, livewrite, has_not_shifted
  use srreheat, only : log_energy_reheat_ingev

  use cndisr, only : cndi_norm_potential
  use cndireheat, only : cndi_x_rreh, cndi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nbeta,nxend

  real(kp) :: alpha,beta,xend,w,bfoldstar,alphamin,alphamax,betamin,betamax,xendmin,xendmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(10) :: alphavalues
  integer, dimension(10) :: nxendvalues

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End
  

  Pstar = powerAmpScalar



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  w=0._kp
  !  w = 1._kp/3._kp

  call delete_file('cndi_predic.dat')
  call delete_file('cndi_nsr.dat')



  npts = 20

!!!!!!!!!!!!!!!!!!!!!!
!!!!    beta=5    !!!!
!!!!!!!!!!!!!!!!!!!!!!

  beta=5._kp

  nalpha=4
  alphavalues(1)=0.01_kp
  alphavalues(2)=0.07_kp
  alphavalues(3)=0.1_kp
  alphavalues(4)=0.12_kp

  nxendvalues(1)=1000
  nxendvalues(2)=200
  nxendvalues(3)=200
  nxendvalues(4)=200

  do j=1,nalpha
     alpha=alphavalues(j)
     nxend=nxendvalues(j)

     xendmin = sqrt(epsilon(1._kp)/2._kp)*(1._kp+beta)/(alpha**2*beta) !to avoid  epsilon_1 < numaccuracy errors
     xendmax = cndi_xendmax(70._kp,alpha,beta)
     xendmin = xendmax/1000._kp


     do k=0,nxend
        xend=xendmin*(xendmax/xendmin)**(real(k,kp)/real(nxend,kp)) !log step


        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = cndi_lnrhoreh_max(alpha,beta,xend,Pstar)

        print *,'alpha=',alpha,'beta=',beta,'xend=',xend,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax


        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = cndi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)

           eps1 = cndi_epsilon_one(xstar,alpha,beta)
           eps2 = cndi_epsilon_two(xstar,alpha,beta)
           eps3 = cndi_epsilon_three(xstar,alpha,beta)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
              cycle
           endif

           
           if ((eps1.lt.1e-5).or.(eps1.gt.0.1) &
                .and.(eps2.lt.0.1).and.(eps2.gt.0.15)) cycle

           call livewrite('cndi_predic.dat',alpha,beta,xend,eps1,eps2,eps3,r,ns,Treh)
          

           call livewrite('cndi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)



        end do
     end do
  end do

!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    beta=0.1    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!

  beta=0.1_kp

  nalpha=6
  alphavalues(1)=0.1_kp
  alphavalues(2)=0.2_kp
  alphavalues(3)=0.3_kp
  alphavalues(4)=0.4_kp
  alphavalues(5)=0.5_kp
  alphavalues(6)=0.6_kp

  nxendvalues(1)=200
  nxendvalues(2)=200
  nxendvalues(3)=100
  nxendvalues(4)=60
  nxendvalues(5)=60
  nxendvalues(6)=60

  do j=1,nalpha
     alpha=alphavalues(j)
     nxend=nxendvalues(j)

     xendmin = sqrt(epsilon(1._kp)/2._kp)*(1._kp+beta)/(alpha**2*beta) !to avoid  epsilon_1 < numaccuracy errors
     xendmax = cndi_xendmax(70._kp,alpha,beta)
     xendmin = xendmax/1000._kp


     do k=0,nxend-1
        xend=xendmin*(xendmax/xendmin)**(real(k,kp)/real(nxend,kp)) !log step


        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = cndi_lnrhoreh_max(alpha,beta,xend,Pstar)

        print *,'alpha=',alpha,'beta=',beta,'xend=',xend,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax


        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = cndi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)

           eps1 = cndi_epsilon_one(xstar,alpha,beta)
           eps2 = cndi_epsilon_two(xstar,alpha,beta)
           eps3 = cndi_epsilon_three(xstar,alpha,beta)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
              cycle
           endif

           
           if ((eps1.lt.1e-5).or.(eps1.gt.0.1) &
                .and.(eps2.lt.0.1).and.(eps2.gt.0.15)) cycle


           call livewrite('cndi_predic.dat',alpha,beta,xend,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('cndi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do
     end do
  end do


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 1e-2
  beta = 5
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = cndi_x_rrad(alpha,beta,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = cndi_epsilon_one(xstar,alpha,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar

     eps1end =  cndi_epsilon_one(xend,alpha,beta)
     VendOverVstar = cndi_norm_potential(xend,alpha,beta)/cndi_norm_potential(xstar,alpha,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = cndi_x_rreh(alpha,beta,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = cndi_x_star(alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program cndimain
