!slow-roll functions for the Logamediate inflation 1 potential ("1"
!means that inflation proceeds from the right to the left)
!
!V(phi) = M^4 x^[4(1-gam)] exp[-beta x^gam]
!
!x = (phi-phi0)/Mp


module lmi1sr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : hypergeom_2F1
  use inftools, only : zbrent
  use lmicommon, only : lmi_alpha
  use lmicommon, only : lmi_norm_potential, lmi_norm_deriv_potential
  use lmicommon, only : lmi_norm_deriv_second_potential
  use lmicommon, only : lmi_epsilon_one_max, lmi_x_potmax
  use lmicommon, only : lmi_epsilon_one, lmi_epsilon_two, lmi_epsilon_three
  use lmicommon, only : lmi_efold_primitive, find_lmi_x_trajectory


  implicit none

  private

  public lmi1_norm_potential
  public lmi1_epsilon_one, lmi1_epsilon_two, lmi1_epsilon_three
  public lmi1_x_endinf, lmi1_efold_primitive, lmi1_x_trajectory
  public lmi1_norm_deriv_potential, lmi1_norm_deriv_second_potential
  public lmi1_betamax
  
contains

!returns V/M^4
  function lmi1_norm_potential(x,gam,beta)    
    implicit none
    real(kp) :: lmi1_norm_potential
    real(kp), intent(in) :: x,gam,beta

    lmi1_norm_potential = lmi_norm_potential(x,gam,beta)

  end function lmi1_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M^4
  function lmi1_norm_deriv_potential(x,gam,beta)
    implicit none
    real(kp) :: lmi1_norm_deriv_potential
    real(kp), intent(in) :: x,gam,beta
    
    lmi1_norm_deriv_potential = lmi_norm_deriv_potential(x,gam,beta)

  end function lmi1_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M^4
  function lmi1_norm_deriv_second_potential(x,gam,beta)
    implicit none
    real(kp) :: lmi1_norm_deriv_second_potential
    real(kp), intent(in) :: x,gam,beta    

    lmi1_norm_deriv_second_potential &
         = lmi_norm_deriv_second_potential(x,gam,beta)
    

  end function lmi1_norm_deriv_second_potential



!epsilon_one(x)
  function lmi1_epsilon_one(x,gam,beta)    
    implicit none
    real(kp) :: lmi1_epsilon_one
    real(kp), intent(in) :: x,gam,beta

    lmi1_epsilon_one = lmi_epsilon_one(x,gam,beta)
    
  end function lmi1_epsilon_one


!epsilon_two(x)
  function lmi1_epsilon_two(x,gam,beta)    
    implicit none
    real(kp) :: lmi1_epsilon_two
    real(kp), intent(in) :: x,gam,beta

    lmi1_epsilon_two = lmi_epsilon_two(x,gam,beta)
    
  end function lmi1_epsilon_two


!epsilon_three(x)
  function lmi1_epsilon_three(x,gam,beta)    
    implicit none
    real(kp) :: lmi1_epsilon_three
    real(kp), intent(in) :: x,gam,beta
   
    lmi1_epsilon_three = lmi_epsilon_three(x,gam,beta)
    
  end function lmi1_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function lmi1_x_endinf(gam,beta)
    implicit none
    real(kp), intent(in) :: gam,beta
    real(kp) :: lmi1_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi1Data

    real(kp) ::xVmax, alpha

    alpha = lmi_alpha(gam)

    xVmax = lmi_x_potmax(gam,beta)

!    maxi = (alpha/(beta*gam))**(1._kp/gam)*(1._kp-epsilon(1._kp))
    maxi = xVmax*(1._kp-epsilon(1._kp))
    mini = epsilon(1._kp)

    lmi1Data%real1 = gam
    lmi1Data%real2 = beta
    
    lmi1_x_endinf = zbrent(find_lmi1_x_endinf,mini,maxi,tolFind,lmi1Data)

  end function lmi1_x_endinf



  function find_lmi1_x_endinf(x,lmi1Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi1Data
    real(kp) :: find_lmi1_x_endinf
    real(kp) :: gam,beta

    gam = lmi1Data%real1
    beta = lmi1Data%real2
    
    find_lmi1_x_endinf = lmi1_epsilon_one(x,gam,beta) - 1._kp
   
  end function find_lmi1_x_endinf

!this is integral[V(phi)/V'(phi) dphi]
  function lmi1_efold_primitive(x,gam,beta)
    implicit none
    real(kp), intent(in) :: x,gam,beta
    real(kp) :: lmi1_efold_primitive

    lmi1_efold_primitive = lmi_efold_primitive(x,gam,beta)
  
  end function lmi1_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function lmi1_x_trajectory(bfold,xend,gam,beta)
    implicit none
    real(kp), intent(in) :: bfold, gam, xend,beta
    real(kp) :: lmi1_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmiData

    real(kp) ::alpha, xVmax

    alpha = lmi_alpha(gam)
    xVmax = lmi_x_potmax(gam,beta)
!    maxi = (alpha/(beta*gam))**(1._kp/gam)*(1._kp-epsilon(1._kp))

    maxi = xVmax*(1._kp-epsilon(1._kp))
    mini = lmi1_x_endinf(gam,beta)*(1._kp+epsilon(1._kp))
  
    lmiData%real1 = gam
    lmiData%real2 = beta
    lmiData%real3 = -bfold + lmi_efold_primitive(xend,gam,beta)
    
    lmi1_x_trajectory = zbrent(find_lmi_x_trajectory,mini,maxi,tolFind,lmiData)
       
  end function lmi1_x_trajectory

! Returns the maximum value of beta such that eps2<eps2max at the top of the potential
  function lmi1_betamax(gam,eps2max)
    implicit none
    real(kp), intent(in) :: gam, eps2max
    real(kp) :: lmi1_betamax

    lmi1_betamax=2._kp**(2-3._kp*gam/2._kp)/ &
                ((1._kp-gam)**(gam/2._kp-1._kp)*gam**(gam/2._kp+1._kp)) &
                *eps2max**(gam/2._kp)

  end function lmi1_betamax
 

end module lmi1sr
