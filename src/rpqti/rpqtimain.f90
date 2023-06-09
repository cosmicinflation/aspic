!test the reheating derivation from slow-roll
program rpqtimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use rpqtisr, only : rpqti_epsilon_one, rpqti_epsilon_two, rpqti_epsilon_three
  use rpqtireheat, only : rpqti_lnrhoreh_max, rpqti_x_star
  use infinout, only : delete_file, livewrite, has_not_shifted
  use srreheat, only : log_energy_reheat_ingev

  use rpqtisr, only : rpqti_norm_potential, rpqti_x_endinf, rpqti_efold_primitive
  use rpqtireheat, only : rpqti_x_rreh, rpqti_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none


  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j,k,l
  integer :: npts,nphi0,nalpha,nbeta

  real(kp) :: phi0,alpha,beta,w,bfoldstar,phi0min,phi0max,alphamin,alphamax,betamin,betamax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) ::x,xmin,xmax

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  integer, parameter :: nvec = 6
  real(kp), dimension(nvec) :: phivec, alphavec
  
  Pstar = powerAmpScalar



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        Calculates the reheating predictions           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call delete_file('rpqti_predic.dat')
  call delete_file('rpqti_nsr.dat')

  call aspicwrite_header('rpqti',labeps12,labnsr,labbfoldreh,(/'beta ','alpha','phi0 '/))
  
  w=0._kp
  !  w = 1._kp/3._kp

  npts = 30

  betamin = -0.9
  betamax = 0.9


  nbeta = 60

  phivec = (/9.0, 9.0, 20.0, 20.0, 50.0, 50.0/)
  alphavec = (/0.1,0.2,-0.2,0.2,-0.2,0.2 /)


  do j=1,nvec
     phi0 = phivec(j)
     alpha = alphavec(j)
     
      do l=0,nbeta
        beta=betamin+(betamax-betamin)*(real(l,kp)/real(nbeta,kp))

        lnRhoRehMin = lnRhoNuc
        xEnd = rpqti_x_endinf(phi0,alpha,beta)
        lnRhoRehMax = rpqti_lnrhoreh_max(phi0,alpha,beta,xend,Pstar)


        print *,'phi0,alpha,beta=',phi0,alpha,beta,'lnRhoRehMin=',lnRhoRehMin, 'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

          lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)


          xstar = rpqti_x_star(phi0,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)


          eps1 = rpqti_epsilon_one(xstar,phi0,alpha,beta)
          eps2 = rpqti_epsilon_two(xstar,phi0,alpha,beta)
          eps3 = rpqti_epsilon_three(xstar,phi0,alpha,beta)


          print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar,'xend=',xend,'xstar=',xstar,'eps1star=',eps1

          logErehGeV = log_energy_reheat_ingev(lnRhoReh)
          Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

          ns = 1._kp - 2._kp*eps1 - eps2
          r =16._kp*eps1

          call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/beta,alpha,phi0/))
          
          if (has_not_shifted(0.002_kp,0.1_kp*log10(eps1),5._kp*eps2)) then
             cycle
          endif

          if ((eps1.lt.1e-5).or.(eps1.gt.0.1) &
               .and.(eps2.lt.0.1).and.(eps2.gt.0.15)) cycle

          call livewrite('rpqti_predic.dat',phi0,alpha,beta,eps1,eps2,eps3,r,ns,Treh)

          call livewrite('rpqti_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        end do

      end do
  
   end do

   call aspicwrite_end()


  write(*,*)
  write(*,*)'Testing Rrad/Rreh'

  lnRradmin=-42
  lnRradmax = 10
  phi0 = 10._kp
  alpha = 0.01_kp
  beta = -0.01_kp

  xend = rpqti_x_endinf(phi0,alpha,beta)
  
  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     xstar = rpqti_x_rrad(phi0,alpha,beta,xend,lnRrad,Pstar,bfoldstar)

     print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

     eps1 = rpqti_epsilon_one(xstar,phi0,alpha,beta)

     !consistency test
     !get lnR from lnRrad and check that it gives the same xstar
     eps1end =  rpqti_epsilon_one(xend,phi0,alpha,beta)
     VendOverVstar = rpqti_norm_potential(xend,phi0,alpha,beta)/rpqti_norm_potential(xstar,phi0,alpha,beta)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     xstar = rpqti_x_rreh(phi0,alpha,beta,xend,lnR,bfoldstar)
     print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar

     !second consistency check
     !get rhoreh for chosen w and check that xstar gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     xstar = rpqti_x_star(phi0,alpha,beta,xend,w,lnRhoReh,Pstar,bfoldstar)
     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar

  enddo


end program rpqtimain
