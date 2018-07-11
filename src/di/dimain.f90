program dimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use infinout, only : delete_file, livewrite

  use dicommon, only : di_direct_x, di_direct_k2, di_k2_nunull, di_k2_potmin
  use dicommon, only : di_norm_parametric_potential, di_norm_deriv_parametric_potential
  use dicommon, only : di_norm_deriv_second_parametric_potential, di_norm_deriv_ln_parametric_potential
  use dicommon, only : di_norm_deriv_third_parametric_potential, di_norm_uplifting
  use dicommon, only : di_parametric_epsilon_one, di_parametric_efold_primitive
  use dicommon, only : di_parametric_epsilon_two, di_parametric_epsilon_three
  use dicommon, only : di_deriv_x, di_deriv_second_x, di_deriv_third_x 
  use displine, only : di_spline_x, di_spline_k2, di_set_splines, di_free_splines
  
  use disr, only : di_x, di_k2, di_norm_potential, di_norm_deriv_potential, di_norm_deriv_second_potential
  use disr, only : di_epsilon_one, di_epsilon_two, di_epsilon_three, di_x_endinf  
  use disr, only : di_k2_epsoneunity, di_k2_trajectory, di_x_trajectory

  use srreheat, only : get_lnrreh_rrad, get_lnrreh_rhow, get_lnrrad_rhow
  use srreheat, only : ln_rho_reheat, ln_rho_endinf, log_energy_reheat_ingev
  use direheat, only : di_k2_star, di_k2_rrad, di_k2_rreh, di_lambda_star, di_lnrhoreh_max
  use direheat, only : di_x_star, di_x_rrad, di_x_rreh

  use infinout, only : aspicwrite_header, aspicwrite_data, aspicwrite_end
  use infinout, only : labeps12, labnsr, labbfoldreh
  
  implicit none

  logical, parameter :: testSpline = .false.
  logical, parameter :: testParametric = .false.
  logical, parameter :: testKahler = .false.




  integer :: i,j,n,npts

  real(kp) :: f, fmin, fmax
  real(kp) :: lambda, ucte

  real(kp) :: k2, lnk2min, lnk2max, k2min, k2max, k2end
  real(kp) :: nu, k2obs, k2star, k2potmin, k2null  
  real(kp) :: xspline, xdirect, dxox, xend, xstar,xpotmin

  real(kp) :: x, xmin, xmax
  real(kp) :: dx, d2x, d3x

  real(kp) :: V, dV, d2V, d3V
  real(kp) :: eps1, eps2, eps3
  real(kp) :: peps1, peps2, peps3

  real(kp) :: w, Pstar, ns, r, Treh, logErehGeV
  real(kp) :: lnRhoReh, lnRhoRehMin, lnRhoRehMax
  real(kp) :: lnRmin, lnRmax, lnR, lnRhoEnd
  real(kp) :: lnRradMin, lnRradMax, lnRrad
  real(kp) :: VendOverVstar, eps1End, bfoldstar




  f=0.01_kp
  print *,'f= ',f
  Pstar = powerAmpScalar
   
  if (testSpline) then
     n=100
     lnk2min = log(1d-8)
     lnk2max = 0._kp

     call di_set_splines()

     do i=2,n
        k2 = exp(lnk2min + real(i-1,kp)*(lnk2max-lnk2min)/real(n-1,kp))
        xdirect = di_direct_x(k2)
        xspline = di_spline_x(k2)
        dxox = 2*(xdirect-xspline)/(xdirect+xspline+epsilon(1._kp))
        print *,'k2= dx/x= ',k2, dxox
     enddo

     call di_free_splines()
  endif
  
  ucte = di_norm_uplifting(f)
  print *,'Uplifting constant=',ucte

  k2null = di_k2_nunull(f)
  print *,'Monopole on at k2= ',k2null

  k2potmin = di_k2_potmin(f)
  print *,'Minimum at k2= ',k2potmin

  xpotmin = di_x(k2potmin)
  print *,'Minimum at x= ',xpotmin

  
  
  if (testParametric) then

     k2min = 1e-3
     k2max = 1-1e-3
     n=100

     call delete_file('parametric_primitive.dat')
     call delete_file('parametric_potential.dat')
     call delete_file('parametric_slowroll.dat')
     call delete_file('parametric_field.dat')

     do i=1,n
        k2 = exp(log(k2min) + real(i-1,kp)*(log(k2max)-log(k2min))/real(n-1,kp))
        
        V = di_norm_parametric_potential(k2,f)
        dV = di_norm_deriv_parametric_potential(k2,f)
        d2V = di_norm_deriv_second_parametric_potential(k2,f)
        d3V = di_norm_deriv_third_parametric_potential(k2,f)

!        print *,'dV/V dlnV= ',dV/V,di_norm_deriv_ln_parametric_potential(k2,f)
        call livewrite('parametric_potential.dat',k2,V,dV,d2V,d3V)

        
        peps1 = di_parametric_epsilon_one(k2,f)
        peps2 = di_parametric_epsilon_two(k2,f)
        peps3 = di_parametric_epsilon_three(k2,f)

        call livewrite('parametric_slowroll.dat',k2,peps1,peps2,peps3)

        x = di_x(k2)
        dx = di_deriv_x(k2)
        d2x = di_deriv_second_x(k2)
        d3x = di_deriv_third_x(k2)

        call livewrite('parametric_field.dat',k2,x,dx,d2x,d3x)

        nu = di_parametric_efold_primitive(k2,f)
        call livewrite('parametric_primitive.dat',k2,nu)

     enddo

  end if



  if (testKahler) then

     call delete_file('potential.dat')
     call delete_file('slowroll.dat')

     lambda = 0.1

     n=1000
     xmin = 0.001
     xmax = 10._kp

     do i=1,n
        x = xmin + real(i-1,kp)*(xmax-xmin)/real(n-1,kp)
        V = di_norm_potential(x,f,lambda)
        dV = di_norm_deriv_potential(x,f,lambda)
        d2V = di_norm_deriv_second_potential(x,f,lambda)

        call livewrite('potential.dat',x,V,dV,d2V)

        eps1 = di_epsilon_one(x,f,lambda)
        eps2 = di_epsilon_two(x,f,lambda)
        eps3 = di_epsilon_three(x,f,lambda)

        call livewrite('slowroll.dat',x,eps1,eps2,eps3)


     enddo

  end if




!  f=1e-3
!  print *,di_lnrhoreh_max(f,Pstar)
!  stop

  call delete_file('di_predic.dat')
  call delete_file('di_nsr.dat')

  call aspicwrite_header('di',labeps12,labnsr,labbfoldreh,(/'f'/))
  
  fmin = 5e-6
  fmax = 0.05
  n = 10

  do j=1,n
     f = exp(log(fmin) + (log(fmax)-log(fmin))*real(j-1,kp)/real(n-1,kp))

     lnRhoRehMin = lnRhoNuc
     lnRhoRehMax = di_lnrhoreh_max(f,Pstar)

     print *,'f= lnRhoRehMin= lnRhoRehMax= ',f,lnRhoRehMin,lnRhoRehMax

     npts = 8
     do i=1,npts

        lnRhoReh = lnRhoRehMin + (lnRhoRehMax-lnRhoRehMin)*real(i-1,kp)/real(npts-1,kp)

        k2star = di_k2_star(f,w,lnRhoReh,Pstar,bfoldstar)
        lambda = di_lambda_star(k2star,f,Pstar)

        print *,'lnRhoReh= ',lnRhoReh, 'lambda= ',lambda, 'bfoldstar= ',bfoldstar
        
        eps1 = di_parametric_epsilon_one(k2star,f)/lambda**2
        eps2 = di_parametric_epsilon_two(k2star,f)/lambda**2
        eps3 = di_parametric_epsilon_three(k2star,f)/lambda**2
       
        logErehGev = log_energy_reheat_ingev(lnRhoReh)
        Treh = 10._kp**( logErehGeV -0.25_kp*log10(acos(-1._kp)**2/30._kp) )

        ns = 1._kp-2._kp*eps1 - eps2
        r = 16._kp*eps1

!for sparing size
        if (abs(ns-1).gt.0.15) cycle
        if (r.lt.1e-10) cycle

        call livewrite('di_predic.dat',f,eps1,eps2,eps3,r,ns,Treh)
        call livewrite('di_nsr.dat',ns,r,abs(bfoldstar),lnRhoReh)

        call aspicwrite_data((/eps1,eps2/),(/ns,r/),(/abs(bfoldstar),lnRhoReh/),(/f/))
        
     enddo
  enddo

  call aspicwrite_end()
  
! Test reheating with lnRrad and lnR

  write(*,*)
  write(*,*)'Testing Rrad/Rreh'
  
  lnRradmin=-40
  lnRradmax = 10

  npts = 3

  do i=1,npts

     lnRrad = lnRradMin + (lnRradMax-lnRradMin)*real(i-1,kp)/real(npts-1,kp)

     k2star = di_k2_rrad(f,lnRrad,Pstar,bfoldstar)
     lambda = di_lambda_star(k2star,f,Pstar)


     print *,'lnRrad=',lnRrad, 'k2star=', k2star, 'lambda= ',lambda, 'bfoldstar= ',bfoldstar

     eps1 = di_parametric_epsilon_one(k2star,f)/lambda**2
     eps2 = di_parametric_epsilon_two(k2star,f)/lambda**2
     eps3 = di_parametric_epsilon_three(k2star,f)/lambda**2     

!consistency test
!get lnR from lnRrad and check that it gives the same xstar
     k2end = di_k2_epsoneunity(f,lambda)
     eps1end =  di_parametric_epsilon_one(k2end,f)/lambda**2
     VendOverVstar = di_norm_parametric_potential(k2end,f) &
          /di_norm_parametric_potential(k2star,f)

     lnRhoEnd = ln_rho_endinf(Pstar,eps1,eps1End,VendOverVstar)
     lnR = get_lnrreh_rrad(lnRrad,lnRhoEnd)
     k2star = di_k2_rreh(f,lnR,Pstar)

     print *,'lnR',lnR,'k2star', k2star

!second consistency check
!get rhoreh for chosen w and check that k2star gotten this way is the same
     w = 0._kp
     lnRhoReh = ln_rho_reheat(w,Pstar,eps1,eps1End,-bfoldstar,VendOverVstar)

     k2star = di_k2_star(f,w,lnRhoReh,Pstar)
     xstar = di_x(k2star)

     print *,'lnR', get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd),'lnRrad' &
          ,get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd),'k2star',k2star



    enddo



end program dimain
