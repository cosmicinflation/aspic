!slow-roll functions for the Logamediate inflation 2 potential ("2"
!means that inflation proceeds from the left to the right)
!
!V(phi) = M^4 x^[4(1-gam)] exp[-beta * x^gam]
!
!x = (phi-phi0)/Mp


module lmi2sr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : hypergeom_2F1
  use inftools, only : zbrent
  use lmicommon, only : lmi_alpha
  use lmicommon, only : lmi_norm_potential, lmi_norm_deriv_potential
  use lmicommon, only : lmi_norm_deriv_second_potential
  use lmicommon, only : lmi_epsilon_one_max, lmi_x_potmax, lmi_x_epsonemax
  use lmicommon, only : lmi_epsilon_one, lmi_epsilon_two, lmi_epsilon_three
  use lmicommon, only : lmi_efold_primitive, find_lmitraj
  implicit none

  private

  public lmi2_norm_potential
  public lmi2_epsilon_one, lmi2_epsilon_two, lmi2_epsilon_three
  public lmi2_efold_primitive, lmi2_x_trajectory
  public lmi2_norm_deriv_potential, lmi2_norm_deriv_second_potential
  public lmi2_xini_min
 
contains

!returns V/M^4
  function lmi2_norm_potential(x,gam,beta)    
    implicit none
    real(kp) :: lmi2_norm_potential
    real(kp), intent(in) :: x,gam,beta

    lmi2_norm_potential = lmi_norm_potential(x,gam,beta)

  end function lmi2_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M^4
  function lmi2_norm_deriv_potential(x,gam,beta)
    implicit none
    real(kp) :: lmi2_norm_deriv_potential
    real(kp), intent(in) :: x,gam,beta
    
    lmi2_norm_deriv_potential = lmi_norm_deriv_potential(x,gam,beta)

  end function lmi2_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M^4
  function lmi2_norm_deriv_second_potential(x,gam,beta)
    implicit none
    real(kp) :: lmi2_norm_deriv_second_potential
    real(kp), intent(in) :: x,gam,beta    

    lmi2_norm_deriv_second_potential &
         = lmi_norm_deriv_second_potential(x,gam,beta)
    

  end function lmi2_norm_deriv_second_potential



!epsilon_one(x)
  function lmi2_epsilon_one(x,gam,beta)    
    implicit none
    real(kp) :: lmi2_epsilon_one
    real(kp), intent(in) :: x,gam,beta

    lmi2_epsilon_one = lmi_epsilon_one(x,gam,beta)
    
  end function lmi2_epsilon_one


!epsilon_two(x)
  function lmi2_epsilon_two(x,gam,beta)    
    implicit none
    real(kp) :: lmi2_epsilon_two
    real(kp), intent(in) :: x,gam,beta

    lmi2_epsilon_two = lmi_epsilon_two(x,gam,beta)
    
  end function lmi2_epsilon_two


!epsilon_three(x)
  function lmi2_epsilon_three(x,gam,beta)    
    implicit none
    real(kp) :: lmi2_epsilon_three
    real(kp), intent(in) :: x,gam,beta
   
    lmi2_epsilon_three = lmi_epsilon_three(x,gam,beta)
    
  end function lmi2_epsilon_three



!returns the minimum value for xin: if eps1<1 in the whole x>xVmax
!interval, returns xVmax, otherwise, returns the highest solution for
!eps1=1
  function lmi2_xini_min(gam,beta)
    implicit none
    real(kp), intent(in) :: gam,beta
    real(kp) :: lmi2_xini_min
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi2Data

    real(kp) ::alpha,xVmax,xeps1max

    alpha = lmi_alpha(gam)

    xVmax = lmi_x_potmax(gam,beta)
    xeps1max = lmi_x_epsonemax(gam,beta)

    if (lmi_epsilon_one_max(gam,beta).lt.1._kp) then
       lmi2_xini_min = xVmax * (1._kp+epsilon(1._kp))

!debug:
       if (.not.(beta.lt.sqrt(2._kp).or.gam.lt.lmi2_gammin(beta))) then
          stop 'lmi2_xini_min: conditional tests WRONG!'
       endif
              
    else

      lmi2_xini_min = lmi2_x_epsoneunity(gam,beta) &
           * (1._kp+epsilon(1._kp))

    endif

  end function lmi2_xini_min



!returns the value xeps1 at which eps1=1 in the domain x>xvmax (if it
!exists, otherwise stop)
  function lmi2_x_epsoneunity(gam,beta)
    implicit none
    real(kp), intent(in) :: gam, beta
    real(kp) :: lmi2_x_epsoneunity
    real(kp) :: alpha, xeps1max
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi2Data

    
    if (lmi_epsilon_one_max(gam,beta).lt.1._kp) then
       stop 'lmi2_x_epsoneunity: no solution for eps1=1'
    endif

    xeps1max = lmi_x_epsonemax(gam,beta)
    alpha = lmi_alpha(gam)

    mini = xeps1max * (1._kp+epsilon(1._kp))
    maxi = 100._kp * max(alpha,(beta*gam)**(1._kp/(1._kp-gam)) &
         ,(alpha*beta*gam)**(1._kp/(2._kp-gam)))
       
    lmi2Data%real1 = gam
    lmi2Data%real2 = beta
    
    lmi2_x_epsoneunity &
         = zbrent(find_lmi2xepsoneunity,mini,maxi,tolFind,lmi2Data)
        
  end function lmi2_x_epsoneunity



  function find_lmi2xepsoneunity(x,lmi2Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi2Data
    real(kp) :: find_lmi2xepsoneunity
    real(kp) :: gam,beta

    gam = lmi2Data%real1
    beta = lmi2Data%real2
    
    find_lmi2xepsoneunity = lmi2_epsilon_one(x,gam,beta) - 1._kp
   
  end function find_lmi2xepsoneunity



!this is integral[V(phi)/V'(phi) dphi]
  function lmi2_efold_primitive(x,gam,beta)
    implicit none
    real(kp), intent(in) :: x,gam,beta
    real(kp) :: lmi2_efold_primitive

    lmi2_efold_primitive = lmi_efold_primitive(x,gam,beta)
    
  end function lmi2_efold_primitive


!this is integral[V(phi)/V'(phi) dphi], approximated in the limit x/x0>>1
  function lmi2_efold_primitive_approximated(x,gam,beta)
    implicit none
    real(kp), intent(in) :: x,gam,beta
    real(kp) :: lmi2_efold_primitive_approximated

    real(kp) ::alpha
    alpha = lmi_alpha(gam)

    if (gam.eq.0._kp) stop 'lmi2_efold_primitive: gam=0!'

    lmi2_efold_primitive_approximated = x**(2._kp-gam) &
         /(beta*gam*(gam-2._kp))

  end function lmi2_efold_primitive_approximated


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function lmi2_x_trajectory(bfold,xend,gam,beta)
    implicit none
    real(kp), intent(in) :: bfold, gam, xend,beta
    real(kp) :: lmi2_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi, xiniMin
    type(transfert) :: lmiData

    real(kp) ::alpha

    alpha = lmi_alpha(gam)
    xiniMin = lmi2_xini_min(gam,beta)

    if (xend.le.xiniMin) then
       write(*,*)'xiniMin= xend= ',xiniMin, xend
       stop 'lmi2_x_trajectory: xend < xiniMin'
    endif

    maxi = xend*(1._kp-epsilon(1._kp))
    mini = xiniMin

    lmiData%real1 = gam
    lmiData%real2 = beta
    lmiData%real3 = -bfold + lmi_efold_primitive(xend,gam,beta)
    
    lmi2_x_trajectory = zbrent(find_lmitraj,mini,maxi,tolFind,lmiData)
       
  end function lmi2_x_trajectory



 
! Returns the minimum value for beta in order to end inflation with
! slow roll violation ( beta>betamin(gam) <=> epsOneMax>1 )
  function lmi2_betamin(gam)
    implicit none
    real(kp), intent(in) :: gam
    real(kp) :: lmi2_betamin
    real(kp) ::alpha

    alpha = lmi_alpha(gam)

    lmi2_betamin=(sqrt(2._kp)*(1._kp-gam)/(alpha*gam))**gam &
         *alpha/(gam*(1._kp-gam))

  end function lmi2_betamin

! Returns the minimum value for gam in order to end inflation with
! slow roll violation ( gam>gammin(beta) <=> epsOneMax>1 )
  function lmi2_gammin(beta)
    implicit none
    real(kp), intent(in) :: beta
    real(kp) :: lmi2_gammin
    real(kp), parameter :: tolFind=tolkp
    type(transfert) :: lmi2Data

    if (beta.lt.sqrt(2._kp)) then
       stop 'lmi2_gammin: beta<sqrt(2): inflation cannot end by slow-roll violation!'
    endif

    lmi2Data%real1 = beta

    lmi2_gammin = zbrent(find_lmi2_gammin,epsilon(1._kp) &
         ,1._kp-epsilon(1._kp),tolFind,lmi2Data)

  end function lmi2_gammin

  function find_lmi2_gammin(x,lmi2Data)
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lmi2Data
    real(kp) :: find_lmi2_gammin
    real(kp)  ::  beta

    beta = lmi2Data%real1
    
    find_lmi2_gammin = lmi2_betamin(x)-beta

  end function find_lmi2_gammin


end module lmi2sr
