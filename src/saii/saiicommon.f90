!common slow-roll functions for string axion inflation I
!
!V(phi) = M^4 [1 - cos(x) + alpha x sin(x)]
!
!x=phi/mu
!
!with no assumptions on alpha and mu
!
module saiicommon
  use infprec, only : kp, tolkp, transfert, pi
  use inftools, only : zbrent, easydverk

  implicit none

  
  private

  public saii_norm_potential, saii_norm_deriv_potential
  public saii_norm_deriv_second_potential
  public saii_epsilon_one, saii_epsilon_two
  public saii_epsilon_three, saii_x_epsoneunity
  public saii_efold_primitive, find_saii_x_trajectory
  public saii_x_potzero, saii_x_derivpotzero, saii_x_potmax

  
contains

  
  function saii_norm_potential(x,alpha,mu)
    implicit none
    real(kp) :: saii_norm_potential
    real(kp), intent(in) :: x,alpha,mu

    saii_norm_potential = 1._kp - cos(x) + alpha*x*sin(x)
    
  end function saii_norm_potential


  
!derivative with respect to x (not phi!)  
  function saii_norm_deriv_potential(x,alpha,mu)
    implicit none
    real(kp) :: saii_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,mu


    saii_norm_deriv_potential = (1._kp + alpha)*sin(x) + alpha*x*cos(x)

  end function saii_norm_deriv_potential

  

!second derivative with respect to x (not phi!)    
  function saii_norm_deriv_second_potential(x,alpha,mu)
    implicit none
    real(kp) :: saii_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,mu
    
    saii_norm_deriv_second_potential = (1._kp+2._kp*alpha)*cos(x) &
         - alpha * x *sin(x)

  end function saii_norm_deriv_second_potential



  
  function saii_epsilon_one(x,alpha,mu)
    implicit none
    real(kp) :: saii_epsilon_one
    real(kp), intent(in) :: x,alpha,mu
    real(kp) :: sinx,cosx
    
    sinx = sin(x)
    cosx = cos(x)
    
    saii_epsilon_one = ( ((1._kp + alpha)*sinx + alpha * x * cosx ) &
         / (1._kp - cosx + alpha * x * sinx) )**2/2._kp/mu/mu
    
  end function saii_epsilon_one
 
  
  
  function saii_epsilon_two(x,alpha,mu)
    implicit none
    real(kp) :: saii_epsilon_two
    real(kp), intent(in) :: x,alpha,mu
    real(kp) :: sinx,cosx
    
    sinx = sin(x)
    cosx = cos(x)
    
    saii_epsilon_two = (2._kp + alpha*(4._kp + alpha + 2._kp*x**2*alpha) &
         - 2._kp*(1._kp + 2._kp*alpha)*cosx - alpha**2*cos(2._kp*x) + 2._kp*x*alpha*sinx) &
         /((1._kp - cosx + x*alpha*sinx)**2)/mu/mu
    
  end function saii_epsilon_two


  
  function saii_epsilon_three(x,alpha,mu)
    implicit none
    real(kp) :: saii_epsilon_three
    real(kp), intent(in) :: x,alpha,mu
    real(kp) :: sinx, cosx, cos2x, sin2x

    cosx = cos(x)
    sinx = sin(x)
    sin2x = sin(2._kp*x)
    cos2x = cos(2._kp*x)
    
    saii_epsilon_three = -(((x*alpha*cosx + (1._kp + alpha)*sinx) &
         * (9._kp*x*alpha**2 - 2._kp*x*alpha*(1._kp + 6._kp*alpha + 2._kp*x**2*alpha**2)*cosx &
         + x*alpha*(2._kp + 3._kp*alpha)*cos2x - 2._kp*sinx - 6._kp*alpha*sinx - 12._kp*alpha**2*sinx &
         - 4._kp*x**2*alpha**2*sinx - 3._kp*alpha**3*sinx + sin2x + 3._kp*alpha*sin2x &
         + 6._kp*alpha**2*sin2x - x**2*alpha**2*sin2x + alpha**3*sin(3*x))) &
         /((1._kp - cosx + x*alpha*sinx)**2*(2._kp + 4._kp*alpha + alpha**2 + 2._kp*x**2*alpha**2 &
         - 2._kp*(1._kp + 2*alpha)*cosx - alpha**2*cos2x + 2._kp*x*alpha*sinx)))/mu/mu

  end function saii_epsilon_three


    
  
!positive field values at which the potential vanishes, the smallest value in ]0,2pi]
  function saii_x_potzero(alpha,mu)
    implicit none
    real(kp) :: saii_x_potzero
    real(kp), intent(in) :: alpha, mu

    real(kp) :: mini, maxi
    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: saiiData

    real(kp), save :: stoalpha = huge(1._kp)
    real(kp), save :: stox = huge(1._kp)
!$omp threadprivate(stoalpha, stox)

    if (alpha.eq.stoalpha) then
       saii_x_potzero = stox
       return
    else
       stoalpha = alpha
    endif
    
    
    if (alpha.gt.0._kp) then

       mini = epsilon(1._kp)
       maxi = 2._kp*pi - epsilon(1._kp)
       saiiData%real1 = alpha
       stox = zbrent(find_saii_x_potzero,mini,maxi,tolFind,saiiData)


    elseif (alpha.ge.-0.5_kp) then
       stox = 2._kp*pi

    else
       mini = epsilon(1._kp)
       maxi = pi - epsilon(1._kp)
       saiiData%real1 = alpha
       stox = zbrent(find_saii_x_potzero,mini,maxi,tolFind,saiiData)

    endif           

    saii_x_potzero = stox
    

  end function saii_x_potzero

  
  function find_saii_x_potzero(x,saiiData)
    implicit none
    real(kp) :: find_saii_x_potzero
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiData
    real(kp) :: alpha

    alpha = saiiData%real1

    find_saii_x_potzero = 1._kp - cos(x) + alpha * x * sin(x)

  end function find_saii_x_potzero
       
 

!non vanishing field values at which the potential is positive and extremal
  function saii_x_derivpotzero(alpha,mu)
    implicit none
    real(kp) :: saii_x_derivpotzero
    real(kp), intent(in) :: alpha,mu

    real(kp) :: mini, maxi
    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: saiiData

    real(kp), save :: stoalpha = huge(1._kp)
    real(kp), save :: stox = huge(1._kp)
!$omp threadprivate(stoalpha,stox)    


    if (alpha.eq.stoalpha) then
       saii_x_derivpotzero = stox
       return
    else
       stoalpha = alpha
    endif
    
    
    if (alpha.gt.0._kp) then

       mini = 0._kp + epsilon(1._kp)
       maxi = pi - epsilon(1._kp)
       saiiData%real1 = alpha
       stox = zbrent(find_saii_x_derivpotzero,mini,maxi,tolFind,saiiData)


    elseif (alpha.eq.0._kp) then
       stox = pi

       
    else
       mini = pi + epsilon(1._kp)
       maxi = 2._kp*pi - epsilon(1._kp)
       saiiData%real1 = alpha       
       stox = zbrent(find_saii_x_derivpotzero,mini,maxi,tolFind,saiiData)

    endif

    saii_x_derivpotzero = stox
    
  end function saii_x_derivpotzero

  
  function find_saii_x_derivpotzero(x,saiiData)
    implicit none
    real(kp) :: find_saii_x_derivpotzero
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiData
    real(kp) :: alpha

    alpha = saiiData%real1

    find_saii_x_derivpotzero = (1._kp + alpha)*sin(x) + alpha*x*cos(x)

  end function find_saii_x_derivpotzero
  


  function saii_x_potmax(alpha,mu)
    implicit none
    real(kp), intent(in) :: alpha,mu
    real(kp) :: saii_x_potmax

    saii_x_potmax = saii_x_derivpotzero(alpha,mu)
    
  end function saii_x_potmax
  



!returns the two roots of epsilon1(x)=1, in the domain for which V>0.
  function saii_x_epsoneunity(alpha,mu)
    implicit none
    real(kp), dimension(2) :: saii_x_epsoneunity
    real(kp), intent(in) :: alpha,mu

    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: saiiData

    real(kp), dimension(2) :: xeps

    real(kp) :: xVzero,xVmax
    real(kp) :: mini, maxi
    integer :: neps2

    saiiData%real1 = alpha
    saiiData%real2 = mu

    xVzero = saii_x_potzero(alpha,mu)
    
    if (alpha.ge.-0.5_kp) then

       xVmax= saii_x_derivpotzero(alpha,mu)

       mini = epsilon(1._kp)
       maxi = xVmax
       saiiData%real3 = 1._kp
       xeps(1) = zbrent(find_saii_x_epsoneunity,mini,maxi,tolFind,saiiData)


       mini = xVmax
       maxi = xVzero - epsilon(1._kp)
       saiiData%real3 = -1._kp
       xeps(2) = zbrent(find_saii_x_epsoneunity,mini,maxi,tolFind,saiiData)

    else

       mini = xVzero + epsilon(1._kp)
       maxi = 2._kp*pi - epsilon(1._kp)
       saiiData%real3 = 1._kp
       xeps(1) = zbrent(find_saii_x_epsoneunity,mini,maxi,tolFind,saiiData)
       
       mini = xVzero + epsilon(1._kp)
       maxi = 2._kp*pi - epsilon(1._kp)
       saiiData%real3 = -1._kp
       xeps(2) = zbrent(find_saii_x_epsoneunity,mini,maxi,tolFind,saiiData)

    end if

    saii_x_epsoneunity = xeps

    
  end function saii_x_epsoneunity

  
  function find_saii_x_epsoneunity(x,saiiData)
    implicit none
    real(kp) :: find_saii_x_epsoneunity
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiData
    real(kp) :: alpha,mu,pm
    real(kp), parameter :: sqr2 = sqrt(2._kp)
    real(kp) :: sinx, cosx, sqrmu

    
    alpha = saiiData%real1
    mu = saiiData%real2
    pm = saiiData%real3

    sinx = sin(x)
    cosx = cos(x)
    sqrmu = sqrt(mu)
    
    find_saii_x_epsoneunity = pm*sqrmu*sqr2*(1._kp - cosx + alpha*x*sinx) &
         - ((1._kp+alpha)*sinx + alpha*x*cosx)/sqrmu
    
  end function find_saii_x_epsoneunity



  
!this is integral[V(phi)/V'(phi) dphi]
  function saii_efold_primitive(x,alpha,mu)
    implicit none
    real(kp), intent(in) :: x, alpha,mu
    real(kp) :: saii_efold_primitive

    type(transfert) :: saiiData
!too long to integrate in QUADPREC
    real(kp), parameter :: tolInt = max(tolkp,epsilon(1._8))
    integer, parameter :: neq = 1

    real(kp) :: xVzero, xpotmax, xvar
    real(kp), dimension(neq) :: yvar

!let us start where the potential vanishes, on the right side compared
!to the maximum
    xVzero = saii_x_potzero(alpha,mu)
    xpotmax = saii_x_derivpotzero(alpha,mu)

    yvar(1) = 0._kp
    
    if (alpha.ge.-0.5_kp) then

       if (x.gt.xpotmax) then
          xvar=xVzero
       else
          xvar=0._kp
       endif

    else

       if (x.gt.xpotmax) then
          xvar=2._kp*pi
       else
          xvar=xVzero
       endif       
    endif


    saiiData%real1 = alpha

    call easydverk(neq,find_saii_efold_primitive,xvar,yvar,x,tolInt,saiiData)

    saii_efold_primitive = yvar(1) * mu * mu
    
  end function saii_efold_primitive

  
  subroutine find_saii_efold_primitive(n,x,y,yprime,saiiData)
    implicit none          
    integer :: n
    real(kp) :: x
    real(kp), dimension(n) :: y, yprime
    type(transfert), optional, intent(inout) :: saiiData
    real(kp) :: alpha,cosx,sinx
!the expansion is accurate up to order 3, and we take a factor of ten
!margin
    real(kp), parameter :: xtaylor = 10._kp*epsilon(1._kp)**(1._kp/3._kp)

    alpha = saiiData%real1

    cosx = cos(x)
    sinx = sin(x)
    
    if (abs(x).lt.xtaylor) then
       if (alpha.eq.-0.5_kp) then
          yprime(1) = 0.25_kp*x
       else
          yprime(1) = 0.5_kp*x
       endif
    else
       yprime(1) = (1._kp - cosx + alpha * x * sinx) &
            / ( (1._kp+alpha)*sinx + alpha * x * cosx )
    endif

  end subroutine find_saii_efold_primitive



  function find_saii_x_trajectory(x,saiiData)
    implicit none
    real(kp) :: find_saii_x_trajectory
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiData

    real(kp) :: alpha, mu, NplusNuend

    alpha = saiiData%real1
    mu = saiiData%real2
    NplusNuend = saiiData%real3

    find_saii_x_trajectory = saii_efold_primitive(x,alpha,mu) - NplusNuend
    
  end function find_saii_x_trajectory

  
end module saiicommon
