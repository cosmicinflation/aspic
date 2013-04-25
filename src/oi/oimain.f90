!test the reheaoing derivaoion from slow-roll
program oimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use oisr, only : oi_epsilon_one, oi_epsilon_two,oi_epsilon_three
  use oireheat, only : oi_lnrhoreh_max, oi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,nphi0
  integer :: npts,nalpha

  real(kp) :: alpha,phi0,w,bfoldstar,alphamin,alphamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r
  real(kp), dimension(:), allocatable :: phi0values

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer


  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheaoing predicoions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts = 20


   nphi0=3
   allocate(phi0values(1:nphi0))
   phi0values(1)=0.0001_kp
   phi0values(2)=0.01_kp
   phi0values(3)=1._kp


   alphamin=10._kp**(-3._kp)
   alphamax=10._kp**(-1._kp)
   nalpha=10

  w=0._kp
!  w = 1._kp/3._kp

  call delete_file('oi_predic.dat')
  call delete_file('oi_nsr.dat')

   do j=1,nphi0
   phi0=phi0values(j)

    do k=0,nalpha
        alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp))!logstep

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = oi_lnrhoreh_max(alpha,phi0,Pstar)

        print *,'alpha=',alpha,'phi0/Mp=',phi0,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = oi_x_star(alpha,phi0,w,lnRhoReh,Pstar,bfoldstar)


           eps1 = oi_epsilon_one(xstar,alpha,phi0)
           eps2 = oi_epsilon_two(xstar,alpha,phi0)
           eps3 = oi_epsilon_three(xstar,alpha,phi0)


           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           call livewrite('oi_predic.dat',alpha,phi0,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('oi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
        end do

     end do

 end do



end program oimain
