!test the reheating derivation from slow-roll
program rpqdimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rpqdisr, only : rpqdi_epsilon_one, rpqdi_epsilon_two, rpqdi_epsilon_three
  use rpqdireheat, only : rpqdi_lnrhoreh_max, rpqdi_x_star
  use infinout, only : delete_file, livewrite, has_not_shifted
  use srreheat, only : log_energy_reheat_ingev

  use rpqdisr, only : rpqdi_norm_potential, rpqdi_x_endinf, rpqdi_efold_primitive
  use rpqdireheat, only : rpqdi_x_rreh, rpqdi_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,l
  integer :: npts,nphi0,nalpha,nbeta

  real(kp) :: phi0,alpha,beta,w,bfoldstar,phi0min,phi0max,alphamin,alphamax,betamin,betamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::x,xmin,xmax, alphaprime, betaprime, normprime
  real(kp) :: alphaprimemin, alphaprimemax, betaprimemax, betaprimemin
  
  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  alphaprimemin = 1e10
  alphaprimemax = -1e10
  betaprimemax = -1e10
  betaprimemin = 1e10
  
  call delete_file('rpqdi_predic.dat')
  call delete_file('rpqdi_nsr.dat')
  call delete_file('rpqdi_prime.dat')

  w=0._kp
  !  w = 1._kp/3._kp

  npts = 5

  alphamin = -0.9
  alphamax = 0.9
  betamin = -0.9
  betamax = 0.9
  phi0min = 5._kp
  phi0max = 1000._kp

  nalpha = 10
  nbeta = 10
  nphi0 = 5

  do j=0,nphi0
     phi0=phi0min*(phi0max/phi0min)**(real(j,kp)/real(nphi0,kp))
    do k=0,nalpha
     alpha=alphamin+(alphamax-alphamin)*(real(k,kp)/real(nalpha,kp))
      do l=0,nbeta
        beta=betamin+(betamax-betamin)*(real(l,kp)/real(nbeta,kp))

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = rpqdi_lnrhoreh_max(phi0,alpha,beta,Pstar)


        print *,'phi0,alpha,beta=',phi0,alpha,beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        call rcpi_conversion(2._kp,alpha,beta,phi0,alphaprime,betaprime,normprime)

        alphaprimemin = min(alphaprime,alphaprimemin)
        alphaprimemax = max(alphaprime, alphaprimemax)
        betaprimemax = max(betaprime, betaprimemax)
        betaprimemin = min(betaprime, betaprimemin)
        
        do i=1,npts

          lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

          xend = rpqdi_x_endinf(phi0,alpha,beta)

          xstar = rpqdi_x_star(phi0,alpha,beta,w,lnRhoReh,Pstar,bfoldstar)


          eps1 = rpqdi_epsilon_one(xstar,phi0,alpha,beta)
          eps2 = rpqdi_epsilon_two(xstar,phi0,alpha,beta)
          eps3 = rpqdi_epsilon_three(xstar,phi0,alpha,beta)


          print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xend=',xend,'xstar=',xstar,'eps1star=',eps1

          logErehGeV = log_energy_reheat_ingev(lnRhoReh)
          Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

          ns = 1._kp - 2._kp*eps1 - eps2
          r =16._kp*eps1

          if (has_not_shifted(0.002_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
             cycle
          endif

          if ((eps1.lt.1e-5).or.(eps1.gt.0.1) &
               .and.(eps2.lt.0.1).and.(eps2.gt.0.15)) cycle

          call livewrite('rpqdi_predic.dat',phi0,alpha,beta,eps1,eps2,eps3,r,ns,Treh)

          call livewrite('rpqdi_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

          call livewrite('rpqdi_prime.dat',alphaprime,betaprime)
          
        end do

      end do
    end do
  end do

  print *,'alphaminmax betaminmax= ',alphaprimemin, alphaprimemax,betaprimemin,betaprimemax


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  phi0 = 10._kp
  alpha = 0.01_kp
  beta = -0.01_kp


  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rpqdi_x_rrad(phi0,alpha,beta,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rpqdi_epsilon_one(xstar,phi0,alpha,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = rpqdi_x_endinf(phi0,alpha,beta)
     eps1end =  rpqdi_epsilon_one(xend,phi0,alpha,beta)
     VendOverVstar = rpqdi_norm_potential(xend,phi0,alpha,beta)/rpqdi_norm_potential(xstar,phi0,alpha,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rpqdi_x_rreh(phi0,alpha,beta,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rpqdi_x_star(phi0,alpha,beta,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

contains

  subroutine rcpi_conversion(p,a,b,mu,alpha,beta,norm)
    implicit none
    real(kp), intent(in) :: p,a,b,mu
    real(kp), intent(out) :: alpha,beta,norm
    real(kp) :: M4om4
    
    M4om4 = 1._kp + 2*(1_kp - a)*log(mu) + 2*(1._kp + b)*log(mu)**2
    
    alpha = -( 2*(1._kp - a) + 4*(1._kp + b)*log(mu))/M4om4
    beta = 2*(1._kp+b)/M4om4
    norm = M4om4/mu**p
    
  end subroutine rcpi_conversion
  
end program rpqdimain
