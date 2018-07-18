!test the reheating derivation from slow-roll
program lmi2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use lmi2sr, only : lmi2_epsilon_one, lmi2_epsilon_two, lmi2_epsilon_three, lmi2_xinimin
  use lmi2reheat, only : lmi2_lnrhoreh_max, lmi2_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  use lmicommon, only : lmi_x_epsonemax, lmi_x_epstwounity,lmi_x_potmax

  use lmi2sr, only : lmi2_norm_potential, lmi2_xendmin
  use lmi2reheat, only : lmi2_x_rreh, lmi2_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,NxEnd
  integer :: npts = 15

  real(kp), dimension(1:6) :: gamValues

  integer, dimension(1:6) :: NxEndValues              

  real(kp) :: xEndMin              !to be specified by lmi2_xinimin
  real(kp) :: xEndMax    

  integer :: Nbeta=10
  real(kp) :: betamin=0.1
  real(kp) :: betamax=10.

  real(kp) :: gam,beta,xEnd,w,bfoldstar,alpha
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End

  Pstar = powerAmpScalar


  call delete_file('lmi2_predic.dat')
  call delete_file('lmi2_nsr.dat')

  call aspicwrite_header('lmi2',labeps12,labnsr,labbfoldreh,(/'xendomin','gamma   ','beta    '/))
  
  !  w = 1._kp/3._kp
  w=0._kp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!    beta=0.1      !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  beta=0.1

  gamValues(1)=0.95_kp
  NxEndValues(1)=15
  gamValues(2)=0.975_kp
  NxEndValues(2)=10
  gamValues(3)=0.995_kp
  NxEndValues(3)=20


  do j=1,3
     gam=gamValues(j)
     NxEnd=nxEndValues(j)

     alpha=4._kp*(1._kp-gam)

     xEndMin=lmi2_xendmin(60._kp,gam,beta)
     xEndMax=xEndMin*20
     
       
     do k=0,NxEnd
        xEnd=xEndMin*(xEndMax/xEndMin)**(real(k,kp)/NxEnd)  !logarithmic step
        ! xEnd=xEndMin+(xEndMax-xEndMin)*(real(k,kp)/NxEnd)  !arithmetic step

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = lmi2_lnrhoreh_max(gam,beta,xEnd,Pstar)

        print *,'gam=',gam,'beta=',beta,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)
           xstar = lmi2_x_star(gam,beta,xEnd,w,lnRhoReh,Pstar,bfoldstar)

           eps1 = lmi2_epsilon_one(xstar,gam,beta)
           eps2 = lmi2_epsilon_two(xstar,gam,beta)
           eps3 = lmi2_epsilon_three(xstar,gam,beta)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)

           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1


           call livewrite('lmi2_predic.dat',gam,beta,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('lmi2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend/xendmin,gam,beta/))

        end do

     end do

  end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!     beta=1       !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  beta=1.

  gamValues(1)=0.6_kp
  NxEndValues(1)=15
  gamValues(2)=0.63_kp
  NxEndValues(2)=10
  gamValues(3)=0.66_kp
  NxEndValues(3)=20

  do j=1,3
     gam=gamValues(j)
     NxEnd=nxEndValues(j)

     alpha=4._kp*(1._kp-gam)
     xEndMin=lmi2_xendmin(60._kp,gam,beta)
     xendmax =  20*xendmin
!     xEndMax=100._kp*max(alpha,(beta*gam)**(1._kp/(1._kp-gam)),(alpha*beta*gam)**(1._kp/(2._kp-gam)))

     do k=0,NxEnd
        xEnd=xEndMin*(xEndMax/xEndMin)**(real(k,kp)/NxEnd)  !logarithmic step
        ! xEnd=xEndMin+(xEndMax-xEndMin)*(real(k,kp)/NxEnd)  !arithmetic step

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = lmi2_lnrhoreh_max(gam,beta,xEnd,Pstar)

        print *,'gam=',gam,'beta=',beta,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)
           xstar = lmi2_x_star(gam,beta,xEnd,w,lnRhoReh,Pstar,bfoldstar)
           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

           eps1 = lmi2_epsilon_one(xstar,gam,beta)
           eps2 = lmi2_epsilon_two(xstar,gam,beta)
           eps3 = lmi2_epsilon_three(xstar,gam,beta)

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)

           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('lmi2_predic.dat',gam,beta,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('lmi2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend/xendmin,gam,beta/))

        end do

     end do

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!     beta=10      !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  beta=10.

  gamValues(1)=0.22_kp
  NxEndValues(1)=10
  gamValues(2)=0.235_kp
  NxEndValues(2)=10
  gamValues(3)=0.24_kp
  NxEndValues(3)=12

  do j=1,3
     gam=gamValues(j)
     NxEnd=nxEndValues(j)

     alpha=4._kp*(1._kp-gam)
     xEndMin=lmi2_xendmin(60._kp,gam,beta)
     xendmax=8*xendmin
!     xEndMax=100._kp*max(alpha,(beta*gam)**(1._kp/(1._kp-gam)),(alpha*beta*gam)**(1._kp/(2._kp-gam)))

     do k=0,NxEnd
        xEnd=xEndMin*(xEndMax/xEndMin)**(real(k,kp)/NxEnd)  !logarithmic step
        ! xEnd=xEndMin+(xEndMax-xEndMin)*(real(k,kp)/NxEnd)  !arithmetic step

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = lmi2_lnrhoreh_max(gam,beta,xEnd,Pstar)

        print *,'gam=',gam,'beta=',beta,'xEnd=',xEnd,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)


           xstar = lmi2_x_star(gam,beta,xEnd,w,lnRhoReh,Pstar,bfoldstar)



           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar


           eps1 = lmi2_epsilon_one(xstar,gam,beta)
           eps2 = lmi2_epsilon_two(xstar,gam,beta)
           eps3 = lmi2_epsilon_three(xstar,gam,beta)


           logErehGeV = log_energy_reheat_ingev(lnRhoReh)

           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('lmi2_predic.dat',gam,beta,xEnd,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('lmi2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/xend/xendmin,gam,beta/))

        end do

     end do

  end do

  call aspicwrite_end()

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  gam = 0.7
  beta = 1.
  xend = lmi2_xendmin(60._kp,gam,beta)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = lmi2_x_rrad(gam,beta,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = lmi2_epsilon_one(xstar,gam,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  lmi2_epsilon_one(xend,gam,beta)
     VendOverVstar = lmi2_norm_potential(xend,gam,beta)/lmi2_norm_potential(xstar,gam,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = lmi2_x_rreh(gam,beta,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = lmi2_x_star(gam,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program lmi2main
