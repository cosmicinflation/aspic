!test the reheating derivation from slow-roll
program sbkimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use sbkisr, only : sbki_epsilon_one, sbki_epsilon_two, sbki_epsilon_three
  use sbkireheat, only : sbki_lnrhoreh_max, sbki_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev

  use sbkisr, only : sbki_norm_potential, sbki_x_endinf
  use sbkisr, only : sbki_efoldmax, sbki_epsilon_one_min
  use sbkisr, only : sbki_xinimax, sbki_x_trajectory, sbki_alphamin, sbki_alphamax
  use sbkireheat, only : sbki_x_rreh, sbki_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  real(kp) :: w,bfoldstar
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp), dimension(1:4) ::alphavalues

  real(kp), parameter :: epsmax = 0.2_kp, efoldNum=80._kp

  real(kp)  :: alpha,alphamin,alphamax,eps1A,eps2A,eps3A,nsA,rA
  real(kp)  :: eps1B,eps2B,eps3B,nsB,rB,xstarA,xstarB
  integer :: Nalpha

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  complex(kp) :: a, b, c, d

  alphavalues(1)=-0.01
  alphavalues(2)=-0.001
  alphavalues(3)=0.001
  alphavalues(4)=0.01

  alphamin=-0.05
  alphamax=0.005
  Nalpha = 200


  Pstar = powerAmpScalar


  call delete_file('sbki_predic.dat')
  call delete_file('sbki_nsr.dat')

!  do j=1,size(alphavalues)
  do j=1,Nalpha

     !alpha=alphavalues(j)
     alpha = alphamin+(alphamax-alphamin)*real(j,kp)/real(Nalpha,kp)

     if ((alpha.ge.sbki_alphamax()).or.(alpha.lt.sbki_alphamin())) cycle

     if (sbki_epsilon_one_min(alpha).gt.epsmax) cycle

     if (sbki_efoldmax(alpha).le.efoldNum) cycle

     w=0._kp

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = sbki_lnrhoreh_max(alpha,Pstar)

     print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = sbki_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = sbki_epsilon_one(xstar,alpha)
        eps2 = sbki_epsilon_two(xstar,alpha)
        eps3 = sbki_epsilon_three(xstar,alpha)


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call livewrite('sbki_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

     end do

  end do


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 0.001
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = sbki_x_rrad(alpha,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = sbki_epsilon_one(xstar,alpha)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     xend = sbki_x_endinf(alpha)
     eps1end =  sbki_epsilon_one(xend,alpha)
     VendOverVstar = sbki_norm_potential(xend,alpha)/sbki_norm_potential(xstar,alpha)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = sbki_x_rreh(alpha,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = sbki_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program sbkimain
