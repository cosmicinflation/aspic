!test the reheating derivation from slow-roll
program nimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use nisr, only : ni_epsilon_one, ni_epsilon_two, ni_epsilon_three
  use nireheat, only : ni_lnrhoreh_max 
  use nireheat, only : ni_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use nisr, only : ni_norm_potential, ni_x_endinf
  use nireheat, only : ni_x_rreh, ni_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat
  
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  real(kp) :: f,w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer


  real(kp), dimension(7) ::fvalues

  real(kp)  ::fmin,fmax,eps1A,eps2A,eps3A,nsA,rA,eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  integer :: nf

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  fvalues(1)=4_kp
  fvalues(2)=5._kp
  fvalues(3)=6._kp
  fvalues(4)=7._kp
  fvalues(5)=10._kp
  fvalues(6)=50._kp
  fvalues(7)=100._kp


  Pstar = powerAmpScalar

  call delete_file('ni_predic.dat')
  call delete_file('ni_nsr.dat')

  call aspicwrite_header('ni',labeps12,labnsr,labbfoldreh,(/'f'/))
  
  do j=1,size(fvalues)

     f=fvalues(j)
     w = 0._kp

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = ni_lnrhoreh_max(f,Pstar)

     print *,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = ni_x_star(f,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,'bfoldstar= ',bfoldstar

        eps1 = ni_epsilon_one(xstar,f)
        eps2 = ni_epsilon_two(xstar,f)
        eps3 = ni_epsilon_three(xstar,f)

        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh =  10._kp**(logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp))


        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('ni_predic.dat',f,eps1,eps2,eps3,r,ns,Treh)

        call livewrite('ni_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/f/))
        
     end do

  end do

  call aspicwrite_end()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Write Data for the summarizing plots !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('ni_predic_summarized.dat') 
  nf=1000
  fmin=10._kp**(0.)
  fmax=10._kp**(2.)
  w=0._kp
  do j=1,nf
     f=fmin*(fmax/fmin)**(real(j,kp)/real(nf,kp))
     lnRhoReh = lnRhoNuc
     xstarA = ni_x_star(f,w,lnRhoReh,Pstar,bfoldstar)
     eps1A = ni_epsilon_one(xstarA,f)
     eps2A = ni_epsilon_two(xstarA,f)
     eps3A = ni_epsilon_three(xstarA,f)
     nsA = 1._kp - 2._kp*eps1A - eps2A
     rA = 16._kp*eps1A
     lnRhoReh = ni_lnrhoreh_max(f,Pstar)
     xstarB = ni_x_star(f,w,lnRhoReh,Pstar,bfoldstar)
     eps1B = ni_epsilon_one(xstarB,f)
     eps2B = ni_epsilon_two(xstarB,f)
     eps3B = ni_epsilon_three(xstarB,f)
     nsB = 1._kp - 2._kp*eps1B - eps2B
     rB =16._kp*eps1B
     call livewrite('ni_predic_summarized.dat',eps1A,eps2A,eps3A,rA,nsA,eps1B,eps2B,eps3B,rB,nsB)
  enddo

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  f = 1.5
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = ni_x_rrad(f,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = ni_epsilon_one(xstar,f)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = ni_x_endinf(f)
     eps1end =  ni_epsilon_one(xend,f)
     VendOverVstar = ni_norm_potential(xend,f)/ni_norm_potential(xstar,f)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = ni_x_rreh(f,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = ni_x_star(f,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo



end program nimain
