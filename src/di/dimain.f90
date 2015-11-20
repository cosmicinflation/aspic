program dimain
  use infprec, only : kp
  use cosmopar, only : lnRhoNuc, powerAmpScalar
  use infinout, only : delete_file, livewrite

  use dicommon, only : di_direct_x, di_direct_k2, di_k2_nunull
  use dicommon, only : di_norm_parametric_potential, di_norm_deriv_parametric_potential
  use dicommon, only : di_norm_deriv_second_parametric_potential
  use dicommon, only : di_norm_deriv_third_parametric_potential
  use dicommon, only :  di_deriv_x, di_deriv_second_x, di_deriv_third_x
  use displine, only : di_spline_x, di_spline_k2, di_set_splines, di_free_splines

  use disr, only : di_x, di_norm_potential, di_norm_deriv_potential, di_norm_deriv_second_potential
  use disr, only : di_epsilon_one, di_epsilon_two, di_epsilon_three, di_x_endinf
  use disr, only : di_parametric_epsilon_one, di_k2_epsoneunity, di_x_trajectory
  use disr, only : di_k2_trajectory, di_x_trajectory, di_parametric_efold_primitive

  implicit none

!test the spline tables giving the correspondance x <--> k2
  logical, parameter :: testSpline = .false.

!test field derivatives
  logical, parameter :: testField = .false.

!dump the potential as a function of field values
  logical, parameter :: testPotential = .true.



  integer :: i,n

  real(kp) :: f, lambda

  real(kp) :: k2, lnk2min, lnk2max, k2min, k2max, k2end
  real(kp) :: nu, k2obs
  real(kp) :: xspline, xdirect, dxox, xend

  real(kp) :: x, xmin, xmax
  real(kp) :: dx, d2x, d3x

  real(kp) :: V, dV, d2V, d3V
  real(kp) :: eps1, eps2, eps3


   f=0.0001_kp
   lambda = 0.01_kp



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


  if (testField) then
     k2min = 0.0001
     k2max = 0.9999
     n=1000

     call delete_file('parametric_field.dat')

     do i=1,n
        k2 = k2min + real(i-1,kp)*(k2max-k2min)/real(n-1,kp)
        x = di_x(k2)
        dx = di_deriv_x(k2)
        d2x = di_deriv_second_x(k2)
        d3x = di_deriv_third_x(k2)

        call livewrite('parametric_field.dat',k2,x,dx,d2x,d3x)

     enddo
  end if


  if (testPotential) then

     k2min = 0.001
     k2max = 0.999
     n=1000

     call delete_file('parametric_potential.dat')
     call delete_file('parametric_trajectory.dat')

     do i=1,n
        k2 = k2min + real(i-1,kp)*(k2max-k2min)/real(n-1,kp)
        
        V = di_norm_parametric_potential(k2,f)
        dV = di_norm_deriv_parametric_potential(k2,f)
        d2V = di_norm_deriv_second_parametric_potential(k2,f)
        d3V = di_norm_deriv_third_parametric_potential(k2,f)

        call livewrite('parametric_potential.dat',k2,V,dV,d2V,d3V)

        nu = di_parametric_efold_primitive(k2,f,lambda)
        call livewrite('parametric_trajectory.dat',k2,nu)


     enddo


     call delete_file('potential.dat')
     call delete_file('slowroll.dat')
    
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

  k2end = di_k2_epsoneunity(f, lambda)
  print *,'k2end= eps1end= ',k2end, di_parametric_epsilon_one(k2end,f,lambda)
  print *,'k2null',di_k2_nunull(f)
  k2obs = di_k2_trajectory(-120._kp,k2end,f,lambda)

  print *,'120 efolds before end: k2=',k2obs,di_x(k2obs)

  

end program dimain
