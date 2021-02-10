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

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20

  real(kp) :: w,bfoldstar,efold
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




  real(kp) :: xmin, xmax, V1,V2,x
  integer :: npt

  
  Pstar = powerAmpScalar



  call delete_file('sbki_potentials.dat')
  call delete_file('sbki_slowrolls.dat')
  


  npt=400

  xmin = 0._kp
  xmax = 30._kp

  do i=1,npt
     x = xmin + real(i-1,kp)*(xmax-xmin)/real(npt-1,kp)        

     V1 = sbki_norm_potential(x,alpha=0.01_kp)
     V2 = sbki_norm_potential(x,alpha=-0.01_kp)

     call livewrite('sbki_potentials.dat',x,V1,V2)

     eps1A = sbki_epsilon_one(x,alpha=0.01_kp)
     eps2A = sbki_epsilon_two(x,alpha=0.01_kp)
     eps3A = sbki_epsilon_three(x,alpha=0.01_kp)

     eps1B = sbki_epsilon_one(x,alpha=-0.01_kp)
     eps2B = sbki_epsilon_two(x,alpha=-0.01_kp)
     eps3B = sbki_epsilon_three(x,alpha=-0.01_kp)
     
     call livewrite('sbki_slowrolls.dat',x,eps1A,eps2A,eps3A,eps1B,eps2B,eps3B)
     

  enddo


  call delete_file('sbki_efoldmax.dat')
  
  npt=300

  alphamin=-0.2
  alphamax=0.1
  
  do i=1,npt
     alpha = alphamin + real(i-1,kp)*(alphamax-alphamin)/real(npt-1,kp) 

     efold = sbki_efoldmax(alpha)

     call livewrite('sbki_efoldmax.dat',alpha,efold)
     
  end do



  
  alphavalues(1)=-0.01
  alphavalues(2)=-0.001
  alphavalues(3)=0.001
  alphavalues(4)=0.01

  alphamin=-0.05
  alphamax=0.005
  Nalpha = 500


  Pstar = powerAmpScalar


  call delete_file('sbki_predic.dat')
  call delete_file('sbki_nsr.dat')

  call aspicwrite_header('sbki',labeps12,labnsr,labbfoldreh,(/'alpha'/))
  
!  do j=1,size(alphavalues)
  do j=1,Nalpha

     !alpha=alphavalues(j)
     alpha = alphamin+(alphamax-alphamin)*real(j,kp)/real(Nalpha,kp)

     if ((alpha.ge.sbki_alphamax()).or.(alpha.lt.sbki_alphamin())) cycle

     if (sbki_epsilon_one_min(alpha).gt.epsmax) cycle

     if (sbki_efoldmax(alpha).le.efoldNum) cycle

     w=0._kp

     lnRhoRehMin = lnRhoNuc
     xEnd = sbki_x_endinf(alpha)
     lnRhoRehMax = sbki_lnrhoreh_max(alpha,xend,Pstar)

     print *,'alpha=',alpha,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        xstar = sbki_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)

        print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xstar=',xstar

        eps1 = sbki_epsilon_one(xstar,alpha)
        eps2 = sbki_epsilon_two(xstar,alpha)
        eps3 = sbki_epsilon_three(xstar,alpha)


        logErehGeV = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp - 2._kp*eps1 - eps2
        r =16._kp*eps1

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/alpha/))
        
        call livewrite('sbki_predic.dat',alpha,eps1,eps2,eps3,r,ns,Treh)

     end do

  end do

  call aspicwrite_end()

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  alpha = 0.001
  xEnd = sbki_x_endinf(alpha)
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = sbki_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = sbki_epsilon_one(xstar,alpha)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  sbki_epsilon_one(xend,alpha)
     VendOverVstar = sbki_norm_potential(xend,alpha)/sbki_norm_potential(xstar,alpha)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = sbki_x_rreh(alpha,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = sbki_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo

end program sbkimain
