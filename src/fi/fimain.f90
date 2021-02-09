!test the reheating derivation from slow-roll
program fimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use fisr, only : fi_epsilon_one, fi_epsilon_two,fi_epsilon_three
  use fireheat, only : fi_lnrhoreh_max, fi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use fisr, only : fi_norm_potential, fi_x_endinf
  use fisr, only : fi_epsilon_one
  use fireheat, only : fi_x_rreh, fi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use fisr, only : fi_efoldmax, fi_epsilon_one_min, fi_x_epsoneunity

  use fisr, only : fi_x_epstwozero

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,ndelta,nn

  real(kp) :: delta,n,w,bfoldstar,deltamin,deltamax,nmin,nmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(:), allocatable ::nvalues

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend
  real(kp), dimension(2) :: xepsone

  real(kp), parameter :: epsmax = 0.2_kp
  real(kp), parameter :: efoldNum = 120._kp

  real(kp) :: Nmax0,Nmax1,Nmax2,Nmax3,Nmax4,Nmax5,Nmax6,Nmax7,Nmax8,Nmax9,Nmax10

  real(kp) :: xmin, xmax, V1, x
  real(kp) :: efold1,efold2,efold3,efold4,efold5,efold6

  integer :: npt

  
  Pstar = powerAmpScalar



  call delete_file('fi_potential.dat')
  call delete_file('fi_slowroll.dat')



  npt=500

  xmin = 0.001_kp
  xmax = 20._kp

  do i=1,npt
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(npt-1,kp)        

     V1 = fi_norm_potential(x,delta=1d-5,n=0._kp)


     call livewrite('fi_potential.dat',x,V1)

     eps1 = fi_epsilon_one(x,delta=1d-5,n=0._kp)
     eps2 = fi_epsilon_two(x,delta=1d-5,n=0._kp)
     eps3 = fi_epsilon_three(x,delta=1d-5,n=0._kp)

     call livewrite('fi_slowroll.dat',x,eps1,eps2,eps3)


  enddo



  
 call delete_file('fi_efoldmax.dat')
  
  npt=150

  deltamin = 1e-10
  deltamax = 0.1
  
  do i=1,npt
     delta = exp(log(deltamin) + real(i-1,kp)*(log(deltamax)-log(deltamin))/real(npt-1,kp))

     efold1 = fi_efoldmax(delta,n=1._kp)
     efold2 = fi_efoldmax(delta,n=2._kp)
     efold3 = fi_efoldmax(delta,n=4._kp)
     efold4 = fi_efoldmax(delta,n=6._kp)
     efold5 = fi_efoldmax(delta,n=8._kp)
     efold6 = fi_efoldmax(delta,n=10._kp)
     
     call livewrite('fi_efoldmax.dat',delta,efold1,efold2,efold3,efold4,efold5,efold6)
     
  end do

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!               Prints Data for Nmax                    !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ndelta=500
deltamin=10._kp**(-12.)
deltamax=10._kp**(-1.)

call delete_file('fi_Nmax.dat')

do k=0,ndelta
  delta=deltamin*(deltamax/deltamin)**(real(k,kp)/real(ndelta,kp))

  Nmax0=0._kp
  Nmax1=0._kp
  Nmax2=0._kp
  Nmax3=0._kp
  Nmax4=0._kp
  Nmax5=0._kp
  Nmax6=0._kp
  Nmax7=0._kp
  Nmax8=0._kp
  Nmax9=0._kp
  Nmax10=0._kp

if (fi_epsilon_one_min(delta,0._kp) .lt. 1._kp) then
  Nmax0=fi_efoldmax(delta,0._kp)
endif

if (fi_epsilon_one_min(delta,1._kp) .lt. 1._kp) then
  Nmax1=fi_efoldmax(delta,1._kp)
endif

if (fi_epsilon_one_min(delta,2._kp) .lt. 1._kp) then
  Nmax2=fi_efoldmax(delta,2._kp)
endif

if (fi_epsilon_one_min(delta,3._kp) .lt. 1._kp) then
  Nmax3=fi_efoldmax(delta,3._kp)
endif

if (fi_epsilon_one_min(delta,4._kp) .lt. 1._kp) then
  Nmax4=fi_efoldmax(delta,4._kp)
endif

if (fi_epsilon_one_min(delta,5._kp) .lt. 1._kp) then
  Nmax5=fi_efoldmax(delta,5._kp)
endif

if (fi_epsilon_one_min(delta,6._kp) .lt. 1._kp) then
  Nmax6=fi_efoldmax(delta,6._kp)
endif

if (fi_epsilon_one_min(delta,7._kp) .lt. 1._kp) then
  Nmax7=fi_efoldmax(delta,7._kp)
endif

if (fi_epsilon_one_min(delta,8._kp) .lt. 1._kp) then
  Nmax8=fi_efoldmax(delta,8._kp)
endif

if (fi_epsilon_one_min(delta,9._kp) .lt. 1._kp) then
  Nmax9=fi_efoldmax(delta,9._kp)
endif

if (fi_epsilon_one_min(delta,10._kp) .lt. 1._kp) then
  Nmax10=fi_efoldmax(delta,10._kp)
endif

  call livewrite('fi_Nmax.dat',delta,Nmax0,Nmax1,Nmax2,Nmax3,Nmax4,Nmax5,Nmax6,Nmax7,Nmax8,Nmax9,Nmax10)

end do



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!
!!      n=0      !!
!!!!!!!!!!!!!!!!!!!

call delete_file('fi_predic_neq0.dat')
call delete_file('fi_nsr_neq0.dat')

call aspicwrite_header('fi',labeps12,labnsr,labbfoldreh,(/'delta','n    '/))

n=0._kp

npts = 10
ndelta = 100

w=0._kp
!  w = 1._kp/3._kp


deltamin=10._kp**(-8.)
deltamax=10._kp**(-3.)

do k=0,ndelta
   delta=deltamin*(deltamax/deltamin)**(real(k,kp)/real(ndelta,kp))

   print*,'delta=',delta,'n=',n

   print*, 'Nmax=',fi_efoldmax(delta,n)
   print*, 'epsilon1min=',fi_epsilon_one_min(delta,n)
   xepsone = fi_x_epsoneunity(delta,n)

   xend = xepsone(1)
   
   print *, 'xepsone', xepsone, fi_epsilon_one(xepsone(1),delta,n),fi_epsilon_one(xepsone(2),delta,n)

    ! Here I have relaxed efoldNum a bit (subtracted 50) in order to produce figures with more spread points, since w=0 does not probe such large DeltaNstar
   if (fi_epsilon_one_min(delta,n).gt.epsmax) cycle
   if (fi_efoldmax(delta,n).lt.(efoldNum-50.)) cycle

   print*, 'hardprior condition passed'


   lnRhoRehMin = lnRhoNuc
   lnRhoRehMax = fi_lnrhoreh_max(delta,n,xend,Pstar)

   print *,lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax


   do i=1,npts

      lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

      xstar = fi_x_star(delta,n,xend,w,lnRhoReh,Pstar,bfoldstar)

      eps1 = fi_epsilon_one(xstar,delta,n)
      eps2 = fi_epsilon_two(xstar,delta,n)
      eps3 = fi_epsilon_three(xstar,delta,n)

      print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1


      logErehGeV = log_energy_reheat_ingev(lnRhoReh)
      Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

      ns = 1._kp - 2._kp*eps1 - eps2
      r =16._kp*eps1

      print*, 'nS=',ns,'r=',r,n


      call livewrite('fi_predic_neq0.dat',delta,n,eps1,eps2,eps3,ns,r,Treh)

      call livewrite('fi_nsr_neq0.dat',ns,r,abs(bfoldstar),lnRhoReh)

      call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/delta,n/))


   end do

end do



!!!!!!!!!!!!!!!!!!!
!!      n=1      !!
!!!!!!!!!!!!!!!!!!!

call delete_file('fi_predic_neq1.dat')
call delete_file('fi_nsr_neq1.dat')

  n=1._kp

  npts = 15
  ndelta=100

  w=0._kp
  !  w = 1._kp/3._kp



  deltamin=10._kp**(-10.)
  deltamax=5*10._kp**(-7.)

  do k=0,ndelta
    delta=deltamin*(deltamax/deltamin)**(real(k,kp)/real(ndelta,kp))

    print*,'delta=',delta,'n=',n

    print*, 'Nmax=',fi_efoldmax(delta,n)
    print*, 'epsilon1min=',fi_epsilon_one_min(delta,n)

    xepsone = fi_x_epsoneunity(delta,n)
    print *, 'xepsone', xepsone, fi_epsilon_one(xepsone(1),delta,n),fi_epsilon_one(xepsone(2),delta,n)

    xend = xepsone(1)
    

    ! Here I have relaxed efoldNum a bit (subtracted 50) in order to produce figures with more spread points, since w=0 does not probe such large DeltaNstar
    if (fi_epsilon_one_min(delta,n).gt.(epsmax*1.)) cycle
    if (fi_efoldmax(delta,n).lt.(efoldNum-50.)) cycle

    print*, 'hardprior condition passed'

    lnRhoRehMin = lnRhoNuc
    lnRhoRehMax = fi_lnrhoreh_max(delta,n,xend,Pstar)

    print *,lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax


    do i=1,npts

      lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

      xstar = fi_x_star(delta,n,xend,w,lnRhoReh,Pstar,bfoldstar)

      eps1 = fi_epsilon_one(xstar,delta,n)
      eps2 = fi_epsilon_two(xstar,delta,n)
      eps3 = fi_epsilon_three(xstar,delta,n)

      print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1


      logErehGeV = log_energy_reheat_ingev(lnRhoReh)
      Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

      ns = 1._kp - 2._kp*eps1 - eps2
      r =16._kp*eps1

      call livewrite('fi_predic_neq1.dat',delta,n,eps1,eps2,eps3,ns,r,Treh)

      call livewrite('fi_nsr_neq1.dat',ns,r,abs(bfoldstar),lnRhoReh)

      call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/delta,n/))
      
    end do

end do


call aspicwrite_end()

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  delta = 10._kp**(-6.)
  n = 0._kp
  xepsone = fi_x_epsoneunity(delta,n)
  xend = xepsone(1)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = fi_x_rrad(delta,n,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = fi_epsilon_one(xstar,delta,n)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  fi_epsilon_one(xend,delta,n)
     VendOverVstar = fi_norm_potential(xend,delta,n)/fi_norm_potential(xstar,delta,n)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = fi_x_rreh(delta,n,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = fi_x_star(delta,n,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo




end program fimain
