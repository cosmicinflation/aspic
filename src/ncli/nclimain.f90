!test the reheating derivation from slow-roll
program nclimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use nclisr, only : ncli_epsilon_one, ncli_epsilon_two, ncli_epsilon_three, ncli_x_endinf
  use nclisr, only : ncli_check_params, ncli_x_epsoneunity, ncli_x_potzero, ncli_x_inflection
  use nclireheat, only : ncli_lnrhoreh_max, ncli_x_star
  use infinout, only : delete_file, livewrite, has_not_shifted
  use srreheat, only : log_energy_reheat_ingev

  use nclisr, only : ncli_norm_potential, ncli_x_epstwozero, ncli_phizeromin
  use nclireheat, only : ncli_x_rreh, ncli_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nphi0,nn

  real(kp) :: alpha,phi0,n,w,bfoldstar,alphamin,alphamax,phi0min,phi0max,nmin,nmax,xEnd
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r,efoldMin, epsMin

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(10) :: alphavalues
  integer, dimension(10) :: nxendvalues

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End
  real(kp) :: xzero, xinf, buf1, buf2
  real(kp), dimension(2) :: xeps2
  real(kp), dimension(3) :: xeps1

  real(kp), dimension(100) :: alphaval
  real(kp) :: xmin,xmax, V1,x

  integer :: npt
  
  Pstar = powerAmpScalar




  npt=250

  xmin = 0._kp
  xmax = 30._kp

  call delete_file('ncli_potential.dat')
  call delete_file('ncli_slowroll.dat')
  
  do i=1,npt
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(npt-1,kp)        

     V1 = ncli_norm_potential(x,alpha=0.001_kp,phi0=10._kp,n=3.0_kp)

     call livewrite('ncli_potential.dat',x,V1)

     eps1 = ncli_epsilon_one(x,alpha=0.001_kp,phi0=10._kp,n=3.0_kp)
     eps2 = ncli_epsilon_two(x,alpha=0.001_kp,phi0=10._kp,n=3.0_kp)
     eps3 = ncli_epsilon_three(x,alpha=0.001_kp,phi0=10._kp,n=3.0_kp)

     call livewrite('ncli_slowroll.dat',x,eps1,eps2,eps3)

  enddo




  

  n=2._kp
  alpha = 0.0001_kp
  phi0 = 0.18_kp



  
  phi0min = ncli_phizeromin(alpha,n)

  print *,'phi0= phi0min= ',phi0,phi0min

  xzero = ncli_x_potzero(alpha,phi0,n)
  xinf = ncli_x_inflection(alpha,phi0,n)
  xeps2 = ncli_x_epstwozero(alpha,phi0,n)
  xeps1 = ncli_x_epsoneunity(alpha,phi0,n)

  print *,'xzero',xzero
  print *,'xinf ',xinf
  print *,'xeps2',xeps2
  print *,'xeps1',xeps1

  do i=1,3
     print *,'eps11',ncli_epsilon_one(xeps1(i),alpha,phi0,n)
  enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  w=0._kp
  !  w = 1._kp/3._kp

  call aspicwrite_header('ncli',labeps12,labnsr,labbfoldreh,(/'phi0 ','alpha','n    '/))
  
  call delete_file('ncli_predic.dat')
  call delete_file('ncli_nsr.dat')

  npts = 15

!!!!!!!!!!!!!!!!!!!!!!
!!!!      n=2     !!!!
!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('ncli_predic_neq2.dat')
  call delete_file('ncli_nsr_neq2.dat')

  n=2._kp


  phi0min = 0.001_kp
  phi0max = 1._kp

  nalpha=2
  alphaval(1)=0.0001_kp
  alphaval(2)=0.000001_kp



  nphi0=100

  efoldMin = 70._kp
  epsMin = 0.1_kp

  do j=1,nalpha
     alpha=alphaval(j)


     do k=0,nphi0
        phi0=phi0min*(phi0max/phi0min)**(real(k,kp)/real(nphi0,kp)) !log step

        if (ncli_check_params(efoldMin,epsMin,alpha,phi0,n)) then

           lnRhoRehMin = lnRhoNuc
           xEnd = ncli_x_endinf(alpha,phi0,n)
           lnRhoRehMax = ncli_lnrhoreh_max(alpha,phi0,n,xend,Pstar)

          print *,'alpha=',alpha,'phi0=',phi0,'n=',n,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

          do i=npts,1,-1
           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)
           xstar = ncli_x_star(alpha,phi0,n,xend,w,lnRhoReh,Pstar,bfoldstar)

           eps1 = ncli_epsilon_one(xstar,alpha,phi0,n)
           eps2 = ncli_epsilon_two(xstar,alpha,phi0,n)
           eps3 = ncli_epsilon_three(xstar,alpha,phi0,n)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )
           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

!           if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
!              cycle
!           endif

!           if ((eps1.lt.1e-10).or.(eps1.gt.0.1) &
!                .and.(eps2.lt.0.1).and.(eps2.gt.0.15)) cycle

           call livewrite('ncli_predic_neq2.dat',alpha,phi0,n,eps1,eps2,ns,r,Treh)
           call livewrite('ncli_nsr_neq2.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/phi0,alpha,n/))

          end do
        end if
     end do
  end do


  n=3._kp

  alphamin = 0.00000001_kp
  alphamax = 0.01_kp
  phi0min = 0.001_kp
  phi0max = 1._kp

  nalpha=2
  nphi0=100

  nalpha=2
  alphaval(1)=0.0001_kp
  alphaval(2)=0.000001_kp
  
  efoldMin = 70._kp
  epsMin = 0.1_kp

  do j=1,nalpha

     alpha = alphaval(j)

     do k=0,nphi0
        phi0=phi0min*(phi0max/phi0min)**(real(k,kp)/real(nphi0,kp)) !log step

        if (ncli_check_params(efoldMin,epsMin,alpha,phi0,n)) then

           lnRhoRehMin = lnRhoNuc
           xEnd = ncli_x_endinf(alpha,phi0,n)
           lnRhoRehMax = ncli_lnrhoreh_max(alpha,phi0,n,xend,Pstar)

          print *,'alpha=',alpha,'phi0=',phi0,'n=',n,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

          do i=npts,1,-1
           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)
           xstar = ncli_x_star(alpha,phi0,n,xend,w,lnRhoReh,Pstar,bfoldstar)

           eps1 = ncli_epsilon_one(xstar,alpha,phi0,n)
           eps2 = ncli_epsilon_two(xstar,alpha,phi0,n)
           eps3 = ncli_epsilon_three(xstar,alpha,phi0,n)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )
           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

!           if (has_not_shifted(0.01_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
!              cycle
!           endif

!           if ((eps1.lt.1e-10).or.(eps1.gt.0.1) &
!                .and.(eps2.lt.0.1).and.(eps2.gt.0.15)) cycle


           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/phi0,alpha,n/))

          end do
        end if
     end do
  end do
  

  call aspicwrite_end()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!      Test Reheating Routines     !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  n=2._kp
  phi0=200._kp
  alpha=0.5_kp
  xEnd = ncli_x_endinf(alpha,phi0,n)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = ncli_x_rrad(alpha,phi0,n,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = ncli_epsilon_one(xstar,alpha,phi0,n)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar

     eps1end =  ncli_epsilon_one(xEnd,alpha,phi0,n)
     VendOverVstar = ncli_norm_potential(xEnd,alpha,phi0,n)/ncli_norm_potential(xstar,alpha,phi0,n)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = ncli_x_rreh(alpha,phi0,n,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = ncli_x_star(alpha,phi0,n,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program nclimain
