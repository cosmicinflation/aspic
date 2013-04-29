!test the reheating derivation from slow-roll
program mlfimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use gmlfisr, only : gmlfi_epsilon_one, gmlfi_epsilon_two, gmlfi_epsilon_three
  use gmlfireheat, only : gmlfi_lnrhoreh_max 
  use gmlfireheat, only : gmlfi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  real(kp) :: p,q,alpha,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  !  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(1:6) ::alphavalues

  real(kp) :: alphamin,alphamax,xstarA,eps1A,eps2A,eps3A,nsA,rA, &
       xstarB,eps1B,eps2B,eps3B,nsB,rB
  integer :: nalpha

  Pstar = powerAmpScalar

  alphavalues(1)=(10._kp)**(-5.)
  alphavalues(2)=(10._kp)**(-3.5)
  alphavalues(3)=(10._kp)**(-3.)
  alphavalues(4)=(10._kp)**(-2.8)
  alphavalues(5)=(10._kp)**(-2.)
  alphavalues(6)=(10._kp)**(0.)

  call delete_file('mlfi_predic.dat')
  call delete_file('mlfi_nsr.dat')

  p = 2._kp 
  q = 2._kp
  w = 0._kp

  do j=1,size(alphavalues)
     alpha=alphavalues(j)

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = gmlfi_lnrhoreh_max(p,q,alpha,Pstar)

     print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = gmlfi_x_star(p,q,alpha,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,'bfoldstar= ',bfoldstar

        eps1 = gmlfi_epsilon_one(xstar,p,q,alpha)
        eps2 = gmlfi_epsilon_two(xstar,p,q,alpha)
        eps3 = gmlfi_epsilon_three(xstar,p,q,alpha)

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh =  10._kp**(logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp))


        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('mlfi_predic.dat',alpha,p,q,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('mlfi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

     end do

  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('mlfi_predic_summarized.dat') 
  nalpha=1000
  alphamin=10._kp**(-5.)
  alphamax=1._kp
  do j=1,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))
     lnRhoReh = lnRhoNuc
     xstarA = gmlfi_x_star(p,q,alpha,w,lnRhoReh,Pstar,bfoldstar)
     eps1A = gmlfi_epsilon_one(xstarA,p,q,alpha)
     eps2A = gmlfi_epsilon_two(xstarA,p,q,alpha)
     eps3A = gmlfi_epsilon_three(xstarA,p,q,alpha)
     nsA = 1._kp - 2._kp*eps1A - eps2A
     rA = 16._kp*eps1A
     lnRhoReh = gmlfi_lnrhoreh_max(p,q,alpha,Pstar)
     xstarB = gmlfi_x_star(p,q,alpha,w,lnRhoReh,Pstar,bfoldstar)
     eps1B = gmlfi_epsilon_one(xstarB,p,q,alpha)
     eps2B = gmlfi_epsilon_two(xstarB,p,q,alpha)
     eps3B = gmlfi_epsilon_three(xstarB,p,q,alpha)
     nsB = 1._kp - 2._kp*eps1B - eps2B
     rB =16._kp*eps1B
     call livewrite('mlfi_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
  enddo


end program mlfimain
