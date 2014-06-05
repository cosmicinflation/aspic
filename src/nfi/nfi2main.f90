program nfi2main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use nfi2sr, only : nfi2_epsilon_one, nfi2_epsilon_two, nfi2_epsilon_three
  use nfi2reheat, only : nfi2_lnrhoreh_max
  use nfi2reheat, only : nfi2_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  use srflow, only : scalar_spectral_index, tensor_to_scalar_ratio

  use nfi2sr, only : nfi2_norm_potential
  use nfi2sr, only : nfi2_amin, nfi2_xendmax, nfi2_numacc_xinimax
  use nfi2sr, only : nfi2_numacc_xendmax, nfi2_numacc_xendmin
  use nfi2reheat, only : nfi2_x_rreh, nfi2_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat
  use nficommon, only : nfi_x_epstwounity,nfi_x_trajectory,nfi_x_epsoneunity
  use srreheat, only : potential_normalization, primscalar

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i
  integer :: npts = 10

  real(kp) :: a,b,w,bfoldstar,y
  real(kp) :: amin,amax,efoldMax,xendmin,xendmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar

  call delete_file('nfi2_predic.dat')
  call delete_file('nfi2_nsr.dat')

  w = 0
  efoldMax=120

  b = 1.8 - 0.4

  do while (b<4)
     
     b = b+0.4

     if (b.lt.2) then
        amin = nfi2_amin(efoldMax,b)
     else
        amin=-huge(1._kp)
     endif

     a=max(amin,-0.5) 

     print *,'a= b= amin= ',a,b,amin
    
     do while (a < -1e-4)

        a = -exp(log(abs(a))-0.5)

        xendmin = nfi2_numacc_xendmin(a,b)

        xendmax = nfi2_numacc_xendmax(efoldMax,a,b)
        
        if (xendmax.lt.xendmin) cycle

        y=1e-2

        do while (y<1)
           
           xend = xendmin + y*(xendmax-xendmin)

           y = exp(log(y)+0.6)

           print *,'a= b= xend= ',a,b,xend
           print *,'xendmin xendmax=',xendmin,xendmax

           lnRhoRehMin = lnRhoNuc

           lnRhoRehMax = nfi2_lnrhoreh_max(a,b,xend,Pstar)

           
           print *,'lnRhoRehMin=',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

           do i=1,npts

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)
              xstar = nfi2_x_star(a,b,xend,w,lnRhoReh,Pstar,bfoldstar)

              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar

              eps1 = nfi2_epsilon_one(xstar,a,b)
              eps2 = nfi2_epsilon_two(xstar,a,b)
              eps3 = nfi2_epsilon_three(xstar,a,b)       

              logErehGeV = log_energy_reheat_ingev(lnRhoReh)
              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              if (abs(ns-1).gt.0.15) cycle
              if (r.lt.1e-10) cycle

              call livewrite('nfi2_predic.dat',a,b,eps1,eps2,eps3,r,ns,Treh)

              call livewrite('nfi2_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           end do           

        enddo       

     enddo

  enddo

! Test reheating with lnRrad and lnR

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'
  
  lnRradmin=-42
  lnRradmax = 10

  b=1.6
  a = -1e-3
  xend = 0.9*nfi2_numacc_xendmax(efoldMax,a,b)

  do i=1,npts

       lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = nfi2_x_rrad(a,b,xend,lnRrad,Pstar,bfoldstar)

       print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

       eps1 = nfi2_epsilon_one(xstar,a,b)
       eps2 = nfi2_epsilon_two(xstar,a,b)
       eps3 = nfi2_epsilon_three(xstar,a,b)
       
!consistency test
!get lnR from lnRrad and check that it gives the same xstar
       eps1end =  nfi2_epsilon_one(xend,a,b)
       VendOverVstar = nfi2_norm_potential(xend,a,b) &
            /nfi2_norm_potential(xstar,a,b)

       lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

       lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
       xstar = nfi2_x_rreh(a,b,xend,lnR,bfoldstar)
       print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar


!second consistency check
!get rhoreh for chosen w and check that xstar gotten this way is the same
       w = 0._kp
       lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)
       
       xstar = nfi2_x_star(a,b,xend,w,lnRhoReh,Pstar,bfoldstar)
       print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
            ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar


                    
    enddo


  end program nfi2main
