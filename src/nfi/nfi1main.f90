program nfi1main
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use nfi1sr, only : nfi1_epsilon_one, nfi1_epsilon_two, nfi1_epsilon_three
  use nfi1reheat, only : nfi1_lnrhoreh_max
  use nfi1reheat, only : nfi1_x_star
  use infinout, only : delete_file, livewrite
  use srreheat, only : log_energy_reheat_ingev
  use srflow, only : scalar_spectral_index, tensor_to_scalar_ratio

  use nfi1sr, only : nfi1_norm_potential, nfi1_x_endinf
  use nfi2sr, only : nfi2_norm_potential
  use nfi1sr, only : nfi1_numacc_amin, nfi1_amax
  use nfi1reheat, only : nfi1_x_rreh, nfi1_x_rrad
  use srreheat, only : get_lnrrad_rreh, get_lnrreh_rrad, ln_rho_endinf
  use srreheat, only : get_lnrrad_rhow, get_lnrreh_rhow, ln_rho_reheat

  use srreheat, only : potential_normalization, primscalar

  
  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none

  
  real(kp) :: Pstar, logErehGeV, Treh

  integer :: i,j
  integer :: npts = 20,n

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

  real(kp) :: x,xmin,xmax,V1,V2

  
  Pstar = powerAmpScalar

  call delete_file('nfi12_potential.dat')
  call delete_file('nfi12_slowroll.dat')

  
  n=250

  xmin = 0.0001_kp
  xmax = 5._kp
  
  do i=1,n
     
     x = exp(log(xmin) + real(i-1,kp)*(log(xmax)-log(xmin))/real(n-1,kp))

     V1 = nfi1_norm_potential(x,a=1._kp,b=1.2_kp)
     V2 = nfi2_norm_potential(x,a=-1._kp,b=1.2_kp)

     call livewrite('nfi12_potential.dat',x,V1,V2)
     
     eps1 = nfi1_epsilon_one(x,a=1._kp,b=1.2_kp)
     eps2 = nfi1_epsilon_two(x,a=1._kp,b=1.2_kp)
     eps3 = nfi1_epsilon_three(x,a=1._kp,b=1.2_kp)

          
     call livewrite('nfi12_slowroll.dat',x,eps1,eps2,eps3)
     
  enddo

  call delete_file('nfi_ab.dat')

  amin = -2._kp
  amax = 2._kp
  
  do i=1,npts
     a = amin+(amax-amin)*real(i-1,kp)/real(npts-1)

     call livewrite('nfi_ab.dat',0.5_kp,1._kp,2._kp,a)

  enddo
     
     
  
  call delete_file('nfi1_predic.dat')
  call delete_file('nfi1_nsr.dat')

  call aspicwrite_header('nfi1',labeps12,labnsr,labbfoldreh,(/'a','b'/))
  
  w = 0
  efoldMax=120


  astep = 0.0001



  bvec=(/3.0,3.5,4.5/)
  
  do j=1,nvec

     b=bvec(j)

     amin = nfi1_numacc_amin(b)
     if (b.lt.2) then
        amax = nfi1_amax(efoldMax,b)
     else
        amax=huge(1._kp)
     endif

     a=max(amin,9e-5)-astep

     print *,'a= b= amin= amax= ',a,b,amin,amax

     do while (a + astep < min(0.1,amax))
 
        a=a+astep
        
        lnRhoRehMin = lnRhoNuc
        xEnd = nfi1_x_endinf(a,b)
        lnRhoRehMax = nfi1_lnrhoreh_max(a,b,xend,Pstar)

        print *,'a= b= ',a,b,'lnRhoRehMin=',lnRhoRehMin,'lnRhoRehMax= ',lnRhoRehMax

        do i=1,npts

           lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

           xstar = nfi1_x_star(a,b,xend,w,lnRhoReh,Pstar,bfoldstar)

           print *,'lnRhoReh',lnRhoReh,' bfoldstar= ',bfoldstar

           eps1 = nfi1_epsilon_one(xstar,a,b)
           eps2 = nfi1_epsilon_two(xstar,a,b)
           eps3 = nfi1_epsilon_three(xstar,a,b)       

           logErehGeV = log_energy_reheat_ingev(lnRhoReh)
           Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

           ns = 1._kp - 2._kp*eps1 - eps2
           r =16._kp*eps1

!           if (abs(ns-1).gt.0.15) cycle
!           if (r.lt.1e-10) cycle

           call livewrite('nfi1_predic.dat',a,b,eps1,eps2,eps3,r,ns,Treh)

           call livewrite('nfi1_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

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

  b=1.8
  a = nfi1_numacc_amin(b)

  xEnd = nfi1_x_endinf(a,b)
  
  do i=1,npts

       lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

       xstar = nfi1_x_rrad(a,b,xend,lnRrad,Pstar,bfoldstar)

       print *,'lnRrad=',lnRrad,' bfoldstar= ',bfoldstar, 'xstar', xstar

       eps1 = nfi1_epsilon_one(xstar,a,b)
       eps2 = nfi1_epsilon_two(xstar,a,b)
       eps3 = nfi1_epsilon_three(xstar,a,b)
       
!consistency test
!get lnR from lnRrad and check that it gives the same xstar
       eps1end =  nfi1_epsilon_one(xend,a,b)
       VendOverVstar = nfi1_norm_potential(xend,a,b) &
            /nfi1_norm_potential(xstar,a,b)

       lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)

       lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
       xstar = nfi1_x_rreh(a,b,xend,lnR,bfoldstar)
       print *,'lnR',lnR, 'bfoldstar= ',bfoldstar, 'xstar', xstar


!second consistency check
!get rhoreh for chosen w and check that xstar gotten this way is the same
       w = 0._kp
       lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)
       
       xstar = nfi1_x_star(a,b,xend,w,lnRhoReh,Pstar,bfoldstar)
       print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
            ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'xstar',xstar


                    
    enddo


  end program nfi1main
