!Slow-roll functions for dual inflation (enjoy!)
!
! V(k2) = M^4{ 1 + Vo(f) - 2(K-E)/(k2 K) - pi/(k2 K K') [nu(k2)]^2 Heavside[nu(k2)] }
! nu(k2)= 1 - 8 sqrt(2)/(pi^2 f) K/sqrt(k2)(E'-K')^2
!
! E=E(k2), K=K(k2), E'=E(1-k2), K'=K(1-k2) are the complete
!elliptic functions, k2=m being the modulus

!Vo(f) is an uplifting term rendering the potential positive and
!depends on the dimensionless parameter f = f0/Lambda. It is defined
!by solving V(k2@minimum) = 0
!
! x = phi/Lambda, the field value in unit of the parameter Lambda,
!stemming from the Kahler potential:
! dx/dk2 = [4 sqrt(2)/pi] sqrt{K K'}/k2^(3/2)
!
!The potential normalization M^4 is completely determined by the model parameters:
! M^4 = f^2 Lambda^4 / pi^2
!

module disr
  use infprec, only : kp,pi,tolkp,transfert
  use inftools, only : zbrent, easydverk
  use dicommon, only : di_norm_parametric_potential, di_norm_deriv_parametric_potential
  use dicommon, only : di_norm_deriv_second_parametric_potential, di_deriv_x
  use dicommon, only : di_norm_deriv_ln_parametric_potential
  use dicommon, only : di_deriv_second_x, di_norm_deriv_third_parametric_potential
  use dicommon, only : di_deriv_third_x, di_k2_nunull, di_k2_potmin
  use dicommon, only : di_parametric_epsilon_one, di_parametric_epsilon_two
  use dicommon, only : di_parametric_epsilon_three, di_parametric_efold_primitive
  implicit none

  private

  logical, parameter :: useKahlerSpline = .true.
  logical, save :: splineSet = .false.

  public di_x, di_k2
  public di_norm_potential, di_epsilon_one, di_epsilon_two, di_epsilon_three
  public di_x_endinf, di_efold_primitive, di_x_trajectory
  public di_norm_deriv_potential, di_norm_deriv_second_potential 
  public di_k2_epsoneunity, di_k2_trajectory


contains

  

!returns the field value x from the parameter k2
  function di_x(k2)
    use dicommon, only : di_direct_x
    use displine, only : di_set_splines, di_spline_x
    use displine, only : k2min, k2max
    implicit none
    real(kp) :: di_x
    real(kp), intent(in) :: k2
!out of the spline, and for k2->0 and k2->1, uses analytical
!integration of first order expansions.
    logical, parameter :: outSplineExp = .true.
    real(kp), parameter :: erfisqrtln4 = 2.3111740399905110619206103738822_kp

    if (.not.useKahlerSpline) then       
       di_x = di_direct_x(k2)
       return
    endif

    if (.not.splineSet) then
       call di_set_splines()
       splineSet = .true.
    endif

    if (outSplineExp) then
!integration of first order expansion of di_deriv_x in k2=1
       if (k2.gt.k2max) then
          di_x = (8._kp*sqrt(2._kp*pi)*erfc(sqrt(log(16._kp/(1._kp - k2)))) &
               -(-1._kp + k2)*sqrt(log(256._kp/(-1._kp + k2)**2)))/sqrt(pi)       
          return
       endif

!integration of first order expansion of di_deriv_x in k2=0
       if (k2.lt.k2min) then
          di_x = (erfisqrtln4 - 8._kp*sqrt(log(4._kp)/pi))/2._kp &
               + (log(65536._kp/k2**4)*log(256._kp/k2**2)**2 &
               - 16._kp*(3._kp + 16._kp*log(2._kp)**2 + log(16._kp) &
               + log(k2)*(-1._kp - 8._kp*log(2._kp) + log(k2)))) &
               / (sqrt(pi*k2)*log(256._kp/k2**2)**2.5_kp)
          return
       end if
    endif

    di_x = di_spline_x(k2)
    
  end function di_x



!returns the parameter k2 from the field value x
  function di_k2(x)
    use dicommon, only : di_direct_k2
    use displine, only : di_set_splines, di_spline_k2
    implicit none
    real(kp) :: di_k2
    real(kp), intent(in) :: x

    if (.not.useKahlerSpline) then       
       di_k2 = di_direct_k2(x)
       return
    endif

    if (.not.splineSet) then
       call di_set_splines()
       splineSet = .true.
    endif

    di_k2 = di_spline_k2(x)

  end function di_k2



!returns M
  function di_potential_normalization(f,lambda)
    real(kp) :: di_potential_normalization
    real(kp), intent(in) :: f, lambda
    real(kp) :: M4

    M4 = f*(lambda**2/pi)**2

    di_potential_normalization = M4**0.25_kp    

  end function di_potential_normalization



!returns Lambda given f and M
  function di_lambda(f,M)
    real(kp) :: di_lambda
    real(kp), intent(in) :: f, M

    di_lambda = M*sqrt(pi/f)
    
  end function di_lambda


!returns V/M^4
  function di_norm_potential(x,f,lambda)
    implicit none
    real(kp) :: di_norm_potential
    real(kp), intent(in) :: x,f,lambda
    real(kp) :: k2    

    k2 = di_k2(x)
   
    di_norm_potential = di_norm_parametric_potential(k2,f)

  end function di_norm_potential



!with respect to x
  function di_norm_deriv_potential(x,f,lambda)
    implicit none
    real(kp) :: di_norm_deriv_potential
    real(kp), intent(in) :: x,f,lambda
    real(kp) :: k2,dx

    if (k2.eq.1._kp) stop 'di_norm_deriv_potential: k2=1, dx=Infinity'

    k2 = di_k2(x)
    dx = di_deriv_x(k2)

    di_norm_deriv_potential = di_norm_deriv_parametric_potential(k2,f)/dx

  end function di_norm_deriv_potential



!with respect to x
  function di_norm_deriv_second_potential(x,f,lambda)
    implicit none
    real(kp) :: di_norm_deriv_second_potential
    real(kp), intent(in) :: x,f,lambda
    real(kp) :: k2, dV, d2V, dx, d2x

    if (k2.eq.1._kp) stop 'di_norm_deriv_second_potential: k2=1, dx=Infinity'

    k2 = di_k2(x)
    dV = di_norm_deriv_parametric_potential(k2,f)
    d2V = di_norm_deriv_second_parametric_potential(k2,f)
    dx = di_deriv_x(k2)
    d2x = di_deriv_second_x(k2)

    di_norm_deriv_second_potential = (d2V - (d2x/dx)*dV)/dx/dx

  end function di_norm_deriv_second_potential


  
  function di_epsilon_one(x,f,lambda)    
    implicit none
    real(kp) :: di_epsilon_one
    real(kp), intent(in) :: x,f,lambda
    real(kp) :: k2

    k2 = di_k2(x)

    di_epsilon_one = di_parametric_epsilon_one(k2,f)/lambda/lambda

  end function di_epsilon_one



  function di_epsilon_two(x,f,lambda)    
    implicit none
    real(kp) :: di_epsilon_two
    real(kp), intent(in) :: x,f,lambda
    real(kp) :: k2

    k2 = di_k2(x)

    di_epsilon_two = di_parametric_epsilon_two(k2,f)/lambda/lambda

  end function di_epsilon_two


  function di_epsilon_three(x,f,lambda)    
    implicit none
    real(kp) :: di_epsilon_three
    real(kp), intent(in) :: x,f,lambda

    real(kp) :: k2

    k2 = di_k2(x)

    di_epsilon_three = di_parametric_epsilon_three(k2,f)/lambda/lambda

  end function di_epsilon_three



  function di_x_endinf(f,lambda)
    implicit none
    real(kp), intent(in) :: f,lambda
    real(kp) :: di_x_endinf

    real(kp) :: k2epsone

    k2epsone = di_k2_epsoneunity(f,lambda)

    di_x_endinf = di_x(k2epsone)
  
  end function di_x_endinf


!returns k2 at which eps1=1
  function di_k2_epsoneunity(f,lambda)
    implicit none
    real(kp) :: di_k2_epsoneunity
    real(kp), intent(in) :: f, lambda

    type(transfert) :: diData
    real(kp), parameter :: tolFind = tolkp
    real(kp) :: mini, maxi , k2potmin

    k2potmin = di_k2_potmin(f)

    mini = tolkp
    maxi = (1._kp-tolkp)*k2potmin
    
    diData%real1 = f
    diData%real2 = lambda

    di_k2_epsoneunity = zbrent(find_di_k2_epsoneunity,mini,maxi,tolFind,diData)   

  end function di_k2_epsoneunity

  
  function find_di_k2_epsoneunity(k2,diData)
    implicit none
    real(kp) :: find_di_k2_epsoneunity
    real(kp), intent(in) :: k2
    type(transfert), optional, intent(inout) :: diData
    real(kp) :: f,lambda

    f = diData%real1
    lambda = diData%real2
    
    find_di_k2_epsoneunity = di_parametric_epsilon_one(k2,f) - lambda**2

  end function find_di_k2_epsoneunity


 
!this is integral[V(phi)/V'(phi) dphi]
  function di_efold_primitive(x,f,lambda)
    implicit none
    real(kp), intent(in) :: x,f,lambda
    real(kp) :: di_efold_primitive

    real(kp) :: k2

    k2 = di_k2(x)

    di_efold_primitive = lambda*lambda*di_parametric_efold_primitive(k2,f)

  end function di_efold_primitive



!returns y at bfold=-efolds before the end of inflation, ie N-Nend
  function di_x_trajectory(bfold,xend,f,lambda)
    implicit none
    real(kp), intent(in) :: bfold, xend, f, lambda
    real(kp) :: di_x_trajectory
    real(kp) :: k2end, k2

    k2end = di_k2(xend)

    k2 = di_k2_trajectory(bfold,k2end,f,lambda)

    di_x_trajectory = di_x(k2)

  end function di_x_trajectory


!the slow-roll trajectory in terms of k2
  function di_k2_trajectory(bfold,k2end,f,lambda)
    implicit none
    real(kp), intent(in) :: bfold, k2end, f, lambda
    real(kp) :: di_k2_trajectory

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: diData


    mini = epsilon(1._kp)
    maxi = k2end


    diData%real1 = f
    diData%real2 = lambda
    diData%real3 = -bfold + lambda*lambda*di_parametric_efold_primitive(k2end,f)

    di_k2_trajectory = zbrent(find_di_k2_trajectory,mini,maxi,tolFind,diData)

  end function di_k2_trajectory


  function find_di_k2_trajectory(k2,diData)
    implicit none
    real(kp), intent(in) :: k2
    type(transfert), optional, intent(inout) :: diData
    real(kp) :: find_di_k2_trajectory
    real(kp) :: f, lambda,NplusNuend

    f = diData%real1
    lambda = diData%real2
    NplusNuend = diData%real3

    find_di_k2_trajectory = lambda*lambda*di_parametric_efold_primitive(k2,f) &
         - NplusNuend 

  end function find_di_k2_trajectory


end module disr
