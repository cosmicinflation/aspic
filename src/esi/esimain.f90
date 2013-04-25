!test the reheating derivation from slow-roll
program esimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use esisr, only : esi_epsilon_one, esi_epsilon_two, esi_epsilon_three
  use esireheat, only : esi_lnrhoreh_max, esi_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use esisr, only : esi_norm_potential, esi_x_endinf
  use esireheat, only : esi_x_rreh, esi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat


  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  real(kp) :: q,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(1:9) ::qvalues

  real(kp)  :: alpha,alphamin,alphamax,eps1A,eps2A,eps3A,nsA,rA
  real(kp)  :: eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  integer :: nalpha

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  qvalues(1)=10._kp**(-3.)
  qvalues(2)=5._kp*10._kp**(-2.)
  qvalues(3)=10._kp**(-1.)
  qvalues(4)=2.5_kp*10._kp**(-1.)
  qvalues(5)=5._kp*10._kp**(-1.)
  qvalues(6)=10._kp**(-1.)
  qvalues(7)=1._kp
  qvalues(8)=1.5_kp
  qvalues(9)=3.5_kp

  Pstar = powerAmpScalar

  call delete_file('esi_predic.dat')
  call delete_file('esi_nsr.dat')

  do j=1,size(qvalues)

     q=qvalues(j)

     w=0._kp

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = esi_lnrhoreh_max(q,Pstar)

     print *,'q=',q,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = esi_x_star(q,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = esi_epsilon_one(xstar,q)
        eps2 = esi_epsilon_two(xstar,q)
        eps3 = esi_epsilon_three(xstar,q)


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('esi_predic.dat',q,w,eps1,eps2,eps3,r,ns,Treh)

     end do
  end do

  do j=1,size(qvalues)

     q=qvalues(j)


     w = -1._kp/3._kp

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = esi_lnrhoreh_max(q,Pstar)

     print *,'q=',q,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = esi_x_star(q,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = esi_epsilon_one(xstar,q)
        eps2 = esi_epsilon_two(xstar,q)
        eps3 = esi_epsilon_three(xstar,q)


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('esi_predic.dat',q,w,eps1,eps2,eps3,r,ns,Treh)

     end do
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('esi_predic_summarized.dat') 
  nalpha=1000
  alphamin=10._kp**(-3.)
  alphamax=10._kp**(1.)
  w=0._kp
  do j=1,nalpha
     alpha=alphamin*(alphamax/alphamin)**(real(j,kp)/real(nalpha,kp))
     lnRhoReh = lnRhoNuc
     xstarA = esi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     eps1A = esi_epsilon_one(xstarA,alpha)
     eps2A = esi_epsilon_two(xstarA,alpha)
     eps3A = esi_epsilon_three(xstarA,alpha)
     nsA = 1._kp - 2._kp*eps1A - eps2A
     rA = 16._kp*eps1A
     lnRhoReh = esi_lnrhoreh_max(alpha,Pstar)
     xstarB = esi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     eps1B = esi_epsilon_one(xstarB,alpha)
     eps2B = esi_epsilon_two(xstarB,alpha)
     eps3B = esi_epsilon_three(xstarB,alpha)
     nsB = 1._kp - 2._kp*eps1B - eps2B
     rB =16._kp*eps1B
     call livewrite('esi_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
  enddo


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  q = 1.
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = esi_x_rrad(q,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = esi_epsilon_one(xstar,q)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = esi_x_endinf(q)
     eps1end =  esi_epsilon_one(xend,q)
     VendOverVstar = esi_norm_potential(xend,q)/esi_norm_potential(xstar,q)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = esi_x_rreh(q,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = esi_x_star(q,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program esimain
