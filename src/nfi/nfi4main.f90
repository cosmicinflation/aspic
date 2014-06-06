program nfi4main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use nfi4sr, only : nfi4_epsilon_one, nfi4_epsilon_two, nfi4_epsilon_three
  use nfi4reheat, only : nfi4_lnrhoreh_max
  use nfi4reheat, only : nfi4_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  use srflow, only : scalar_spectral_index, tensor_to_scalar_ratio

  use nfi4sr, only : nfi4_norm_potential
  use nfi4sr, only : nfi4_xendmin, nfi4_numacc_xendmax
  use nfi4reheat, only : nfi4_x_rreh, nfi4_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat
  use srreheat, only : potential_normalization, primscalar

  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i
  integer :: npts = 10

  real(kp) :: a,b,w,bfoldstar,y
  real(kp) :: astep, bstep, ystep
  real(kp) :: amin,amax,efoldMax,xendmin,xendmax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  Pstar = powerAmpScalar

  call delete_file('nfi4_predic.dat')
  call delete_file('nfi4_nsr.dat')

  w = 0
  efoldMax=120

!region a>0 and 0<b<1
  bstep = 0.2
  astep = 0.3
  ystep = 0.2


  b = 0.001 - bstep

  do while (b+bstep<1)
     
     b = b+bstep
         
     a = 0.05

     do while (a + astep < 2)

        a = a+astep

        xendmin = nfi4_xendmin(efoldMax,a,b)
        
        xendmax = nfi4_numacc_xendmax(a,b)

        if (xendmax.le.xendmin) cycle

        y=1e-2

        do while (y<1)
           
           xend = xendmin + y*(xendmax-xendmin)

           y = exp(log(y)+ystep)

           print *,'a= b= xend= ',a,b,xend
           print *,'xendmin xendmax=',xendmin,xendmax

           lnRhoRehMin = lnRhoNuc

           lnRhoRehMax = nfi4_lnrhoreh_max(a,b,xend,Pstar)

           
           print *,'lnRhoRehMin=',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

           do i=1,npts

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)
              xstar = nfi4_x_star(a,b,xend,w,lnRhoReh,Pstar,bfoldstar)

              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar

              eps1 = nfi4_epsilon_one(xstar,a,b)
              eps2 = nfi4_epsilon_two(xstar,a,b)
              eps3 = nfi4_epsilon_three(xstar,a,b)       

              logErehGeV = log_energy_reheat_ingev(lnRhoReh)
              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              if (abs(ns-1).gt.0.15) cycle
              if (r.lt.1e-10) cycle

              call livewrite('nfi4_predic.dat',a,b,eps1,eps2,eps3,r,ns,Treh)

              call livewrite('nfi4_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           end do           

        enddo       

     enddo

  enddo


!region a<0 and b<0
  bstep = 0.2
  astep = 0.3
  ystep = 0.2


  b = -2 - bstep

  do while (b+bstep<0)
     
     b = b+bstep
         
     a = -3

     do while (a + astep < 0)

        a = a+astep

        xendmin = nfi4_xendmin(efoldMax,a,b)
        
        xendmax = nfi4_numacc_xendmax(a,b)

        print *,'test',xendmin,xendmax

        if (xendmax.le.xendmin) cycle

        y=1e-2

        do while (y<1)
           
           xend = xendmin + y*(xendmax-xendmin)

           y = exp(log(y)+ystep)

           print *,'a= b= xend= ',a,b,xend
           print *,'xendmin xendmax=',xendmin,xendmax

           lnRhoRehMin = lnRhoNuc

           lnRhoRehMax = nfi4_lnrhoreh_max(a,b,xend,Pstar)

           
           print *,'lnRhoRehMin=',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

           do i=1,npts

              lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)
              xstar = nfi4_x_star(a,b,xend,w,lnRhoReh,Pstar,bfoldstar)

              print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar

              eps1 = nfi4_epsilon_one(xstar,a,b)
              eps2 = nfi4_epsilon_two(xstar,a,b)
              eps3 = nfi4_epsilon_three(xstar,a,b)       

              logErehGeV = log_energy_reheat_ingev(lnRhoReh)
              Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

              ns = 1._kp - 2._kp*eps1 - eps2
              r =16._kp*eps1

              if (abs(ns-1).gt.0.15) cycle
              if (r.lt.1e-10) cycle

              call livewrite('nfi4_predic.dat',a,b,eps1,eps2,eps3,r,ns,Treh)

              call livewrite('nfi4_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           end do           

        enddo       

     enddo

  enddo




! Test reheating with lnRrad and lnR

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'
  
  lnRradmin=-42
  lnRradmax = 10

  b=0.8
  a = 0.2

  xend = 2*nfi4_xendmin(efoldMax,a,b)

  do i=1,npts

       lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = nfi4_x_rrad(a,b,xend,lnRrad,Pstar,bfoldstar)

       print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

       eps1 = nfi4_epsilon_one(xstar,a,b)
       eps2 = nfi4_epsilon_two(xstar,a,b)
       eps3 = nfi4_epsilon_three(xstar,a,b)
       
!consistency test
!get lnR from lnRrad and check that it gives the same xstar
       eps1end =  nfi4_epsilon_one(xend,a,b)
       VendOverVstar = nfi4_norm_potential(xend,a,b) &
            /nfi4_norm_potential(xstar,a,b)

       lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

       lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
       xstar = nfi4_x_rreh(a,b,xend,lnR,bfoldstar)
       print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar


!second consistency check
!get rhoreh for chosen w and check that xstar gotten this way is the same
       w = 0._kp
       lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)
       
       xstar = nfi4_x_star(a,b,xend,w,lnRhoReh,Pstar,bfoldstar)
       print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
            ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar


                    
    enddo


  end program nfi4main
