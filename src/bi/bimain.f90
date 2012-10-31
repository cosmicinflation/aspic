!test the reheating derivation from slow-roll
program kklttest  
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use kkltreheat, only : kklt_x_star, kklt_lnrhoend
  use kkltsrevol, only : kklt_x_endinf, kklt_epsilon_one, kklt_epsilon_two
  use infinout, only : delete_file, livewrite
  implicit none
 
  real(kp) :: Pstar

  integer :: i
  integer :: npts = 20

  real(kp) :: mu,p,wreh,bfold
  real(kp) :: lnRhoReh,chi,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(3) :: vecbuffer

  logical, parameter :: display = .true.

  p = 4._kp
  mu = 1._kp
  wreh = -0.3_kp

  Pstar = powerAmpScalar


  call delete_file('kkltreh_eps2eps1.dat')
  call delete_file('kkltreh_nsr.dat')

!chi stands for phi/mu
 
  lnRhoRehMin = lnRhoNuc
  lnRhoRehMax = kklt_lnrhoend(p,mu,Pstar)

  print *,'lnRhoRehMin= lnRhoRehMax= ',lnRhoRehMin,lnRhoRehMax

  do i=1,npts

     lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

     chi = kklt_x_star(p,mu,wreh,lnRhoReh,Pstar,bfold)     
     eps1 = kklt_epsilon_one(chi,p,mu)
     eps2 = kklt_epsilon_two(chi,p,mu)

     if (display) print *,'lnRhoReh= N*= ',lnRhoReh,abs(bfold),eps1

     call livewrite('kkltreh_eps2eps1.dat',eps2,eps1,abs(bfold),lnRhoReh)
   
     ns = 1._kp - 2._kp*eps1 - eps2
     r =16._kp*eps1

     call livewrite('kkltreh_nsr.dat',ns,r,abs(bfold),lnRhoReh)

  end do

 
end program kklttest
