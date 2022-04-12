!test the reheating derivation from slow-roll
program saaimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use saaisr, only : saai_epsilon_one, saai_epsilon_two, saai_epsilon_three
  use saaireheat, only : saai_lnrhoreh_max, saai_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use saaisr, only : saai_norm_potential, saai_x_endinf, saai_x_trajectory
  use saaireheat, only : saai_x_rreh, saai_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,n
  integer :: npts = 20

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp)  :: alpha,alphamin,alphamax,eps1A,eps2A,eps3A,nsA,rA
  real(kp)  :: eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  integer :: Nalpha

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  real(kp) :: x, xmin, xmax, V1

  Pstar = powerAmpScalar

  call delete_file('saai_potential.dat')
  call delete_file('saai_slowroll.dat')

  n=100

  xmin = 1e-2_kp
  xmax = 5._kp
  
  do i=1,n
     x = exp(log(xmin) + real(i-1,kp)*(log(xmax)-log(xmin))/real(n-1,kp))

     V1 = saai_norm_potential(x,alpha=0.1_kp)


     call livewrite('saai_potential.dat',x,V1)

     if (x.eq.0._kp) cycle
     
     eps1 = saai_epsilon_one(x,alpha=0.1_kp)
     eps2 = saai_epsilon_two(x,alpha=0.1_kp)
     eps3 = saai_epsilon_three(x,alpha=0.1_kp)
          
     call livewrite('saai_slowroll.dat',x,eps1,eps2,eps3)
     
  enddo
  



  
  alphamin=0.001
  alphamax=100.
  Nalpha = 1000




  call delete_file('saai_predic.dat')
  call delete_file('saai_nsr.dat')

  call aspicwrite_header('saai',labeps12,labnsr,labbfoldreh,(/'alpha'/))
  
!  do j=1,size(alphavalues)
  do j=1,Nalpha

     !alpha=alphavalues(j)
     alpha = exp(log(alphamin)+(log(alphamax)-log(alphamin))*real(j,kp)/real(Nalpha,kp))

     w=0._kp

     lnRhoRehMin = lnRhoNuc
     xEnd = saai_x_endinf(alpha)
     lnRhoRehMax = saai_lnrhoreh_max(alpha,xend,Pstar)

     print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = saai_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = saai_epsilon_one(xstar,alpha)
        eps2 = saai_epsilon_two(xstar,alpha)
        eps3 = saai_epsilon_three(xstar,alpha)


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha/))
        
        call livewrite('saai_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

     end do

  end do

  call aspicwrite_end()

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 2._kp

  xEnd = saai_x_endinf(alpha)

  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = saai_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = saai_epsilon_one(xstar,alpha)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  saai_epsilon_one(xend,alpha)
     VendOverVstar = saai_norm_potential(xend,alpha)/saai_norm_potential(xstar,alpha)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = saai_x_rreh(alpha,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = saai_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program saaimain
