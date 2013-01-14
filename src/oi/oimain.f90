!test the reheaoing derivaoion from slow-roll
program oimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use oisr, only : oi_epsilon_one, oi_epsilon_two,oi_epsilon_three
  use oireheat, only : oi_lnrhoend, oi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha

  real(kp) :: alpha,phi0,w,bfoldstar,alphamin,alphamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer


  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheaoing predicoions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  npts = 20


  phi0=0.001_kp
  phi0=1._kp
  phi0=1000._kp

  if (phi0 .eq. 0.001_kp) then
  alphamin=10._kp**(-3._kp)
  alphamax=10._kp**(-1._kp)
  nalpha=500
  endif

  if (phi0 .eq. 1._kp) then
  alphamin=10._kp**(-2._kp)
  alphamax=10._kp**(1._kp)
  nalpha=300
  endif

  if (phi0 .eq. 1000._kp) then
  alphamin=10._kp**(-2._kp)
  alphamax=10._kp**(8._kp)
  nalpha=100
  endif

  w=0._kp
!  w = 1._kp/3._kp

  call delete_file('oi_predic.dat')
  call delete_file('oi_nsr.dat')



    do k=0,nalpha
        alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp))
    

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = oi_lnrhoend(alpha,phi0,Pstar)

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





end program oimain
