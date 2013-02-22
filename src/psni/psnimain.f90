!test the reheating derivation from slow-roll
program psnimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use psnisr, only : psni_epsilon_one, psni_epsilon_two, psni_epsilon_three
  use psnireheat, only : psni_lnrhoend, psni_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev


  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k
  integer :: npts,nalpha,nf

  real(kp) :: alpha,f,w,bfoldstar,alphamin,alphamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r
  real(kp), dimension(:), allocatable :: fvalues


  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer


  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call delete_file('psni_predic.dat')
  call delete_file('psni_nsr.dat')
  
  npts = 20 
  nalpha=200
  w=0._kp
!  w = 1._kp/3._kp

   
   nf=3
   allocate(fvalues(1:3))
   fvalues(3)=0.001_kp
   fvalues(2)=0.1_kp
   fvalues(1)=10._kp
   
   do j=1,nf
   f=fvalues(j)

     alphamax=10._kp**(-1._kp)*f**2
     alphamin=10._kp**(-8._kp)*f**2

    do k=1,nalpha
       alpha=alphamin*(alphamax/alphamin)**(real(k,kp)/real(nalpha,kp)) !log step
     
      lnRhoRehMin = lnRhoNuc
      lnRhoRehMax = psni_lnrhoend(alpha,f,Pstar)

      print *,'alpha=',alpha,'f=',f,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

      do i=1,npts

       lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = psni_x_star(alpha,f,w,lnRhoReh,Pstar,bfoldstar)


       eps1 = psni_epsilon_one(xstar,alpha,f)
       eps2 = psni_epsilon_two(xstar,alpha,f)
       eps3 = psni_epsilon_three(xstar,alpha,f)

       print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar,'eps1star=',eps1

       logErehGeV = log_energy_reheat_ingev(lnRhoReh)
       Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

       ns = 1._kp - 2._kp*eps1 - eps2
       r =16._kp*eps1

       call livewrite('psni_predic.dat',alpha,f,eps1,eps2,eps3,r,ns,Treh)

       call livewrite('psni_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)
  
    end do

  end do

 end do



end program psnimain
