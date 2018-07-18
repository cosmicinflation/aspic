program nfi3main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use nfi3sr, only : nfi3_epsilon_one, nfi3_epsilon_two, nfi3_epsilon_three
  use nfi3reheat, only : nfi3_lnrhoreh_max
  use nfi3reheat, only : nfi3_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  use srflow, only : scalar_spectral_index, tensor_to_scalar_ratio

  use nfi3sr, only : nfi3_norm_potential, nfi3_x_endinf, nfi3_numacc_xinimax
  use nfi3sr, only : nfi3_numacc_absamax
  use nfi3reheat, only : nfi3_x_rreh, nfi3_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat
  use nficommon, only : nfi_numacc_x_potbig
  use srreheat, only : potential_normalization, primscalar

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 15

  real(kp) :: a,b,w,bfoldstar
  real(kp) :: astep, bstep
  real(kp) :: amin,amax,efoldMax
  real(kp) :: lnRhoReh,xstar,eps1,eps2,eps3,ns,r

  real(kp) :: lnRhoRehMin, lnRhoRehMax
  real(kp), dimension(2) :: vecbuffer

  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, xend

  integer, parameter :: nvec = 3
  real(kp), dimension(nvec) :: bvec

  
  Pstar = powerAmpScalar

  call delete_file('nfi3_predic.dat')
  call delete_file('nfi3_nsr.dat')

  
  call aspicwrite_header('nfi3',labeps12,labnsr,labbfoldreh,(/'a','b'/))
  
  w = 0
  efoldMax=120

  astep = 0.03

!region a<0 and 0<b<1

  bvec=(/0.05, 0.25, 0.30/)

  do j=1,nvec

     b=bvec(j)

     amin = -nfi3_numacc_absamax(b)
     
     a = max(-10.,amin) - astep

     print *,'a= b= amin= ',a,b,amin

    

     do while (a+astep < -1e-4)
 
        a=a+astep

        print *,'xend',nfi3_x_endinf(a,b)
        print *,'xinimax',nfi3_numacc_xinimax(a,b)
        print *,'xpotbig',nfi_numacc_x_potbig(a,b)

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = nfi3_lnrhoreh_max(a,b,Pstar)

        print *,'a= b= ',a,b,'lnRhoRehMin=',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = nfi3_x_star(a,b,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar

           eps1 = nfi3_epsilon_one(xstar,a,b)
           eps2 = nfi3_epsilon_two(xstar,a,b)
           eps3 = nfi3_epsilon_three(xstar,a,b)       

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           if (abs(ns-1).gt.0.15) cycle
           if (r.lt.1e-10) cycle

           call livewrite('nfi3_predic.dat',a,b,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('nfi3_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/a,b/))
           
        end do

     enddo
  enddo



!region a>0 and b<0

  astep = 0.03

  bvec=(/-0.1,-1.1,-2.1/)

  do j=1,nvec

     b=bvec(j)

     amax = nfi3_numacc_absamax(b)

     a = min(amax,1e-4) - astep

     print *,'a= b=  amax= ',a,b,amax



     do while (a + astep < min(4._kp,amax))

        a=a+astep

        print *,'xend',nfi3_x_endinf(a,b)
        print *,'xinimax',nfi3_numacc_xinimax(a,b)
        print *,'xpotbig',nfi_numacc_x_potbig(a,b)

        lnRhoRehMin = lnRhoNuc
        lnRhoRehMax = nfi3_lnrhoreh_max(a,b,Pstar)

        print *,'a= b= ',a,b,'lnRhoRehMin=',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = nfi3_x_star(a,b,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar

           eps1 = nfi3_epsilon_one(xstar,a,b)
           eps2 = nfi3_epsilon_two(xstar,a,b)
           eps3 = nfi3_epsilon_three(xstar,a,b)       

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

           if (abs(ns-1).gt.0.15) cycle
           if (r.lt.1e-10) cycle

           call livewrite('nfi3_predic.dat',a,b,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('nfi3_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

           call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/a,b/))
           
        end do

     enddo
  enddo


  call aspicwrite_end()


! Test reheating with lnRrad and lnR

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'
  
  lnRradmin=-42
  lnRradmax = 10

  b=-4
  a = 2


  do i=1,npts

       lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = nfi3_x_rrad(a,b,lnRrad,Pstar,bfoldstar)

       print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

       eps1 = nfi3_epsilon_one(xstar,a,b)
       eps2 = nfi3_epsilon_two(xstar,a,b)
       eps3 = nfi3_epsilon_three(xstar,a,b)
       
!consistency test
!get lnR from lnRrad and check that it gives the same xstar
       xend = nfi3_x_endinf(a,b)
       eps1end =  nfi3_epsilon_one(xend,a,b)
       VendOverVstar = nfi3_norm_potential(xend,a,b) &
            /nfi3_norm_potential(xstar,a,b)

       lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

       lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
       xstar = nfi3_x_rreh(a,b,lnR,bfoldstar)
       print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar


!second consistency check
!get rhoreh for chosen w and check that xstar gotten this way is the same
       w = 0._kp
       lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)
       
       xstar = nfi3_x_star(a,b,w,lnRhoReh,Pstar,bfoldstar)
       print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
            ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar


                    
    enddo


  end program nfi3main
