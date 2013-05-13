!slow-roll functions for the rpower law inflation potential
!
!V(phi) = M^4 exp(- alpha x )
!
!x = phi/Mp

module plisr
  use infprec, only : kp
  implicit none

  private

  real(kp), parameter :: junkXend = 100._kp

  public  pli_norm_potential, pli_epsilon_one, pli_epsilon_two, pli_epsilon_three
  public  pli_efold_primitive, pli_x_trajectory, pli_x_endinf
  public  pli_norm_deriv_potential, pli_norm_deriv_second_potential
 
contains
!returns V/M^4
  function pli_norm_potential(x,alpha)
    implicit none
    real(kp) :: pli_norm_potential
    real(kp), intent(in) :: x,alpha

    pli_norm_potential = exp(-alpha*x)
  end function pli_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function pli_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: pli_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

   pli_norm_deriv_potential = -alpha*exp(-alpha*x)

  end function pli_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function pli_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: pli_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    pli_norm_deriv_second_potential = alpha**2*exp(-alpha*x)

  end function pli_norm_deriv_second_potential



!epsilon_one(x)
  function pli_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: pli_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    pli_epsilon_one = alpha**2/2._kp
    
  end function pli_epsilon_one


!epsilon_two(x)
  function pli_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: pli_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    pli_epsilon_two = 0._kp
    
  end function pli_epsilon_two


!epsilon_three(x)
  function pli_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: pli_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    pli_epsilon_three = 0._kp
    
  end function pli_epsilon_three


  function pli_x_endinf(alpha)
    implicit none
    real(kp) :: pli_x_endinf
    real(kp), intent(in) :: alpha

    pli_x_endinf = junkXend

  end function pli_x_endinf



!this is integral[V(phi)/V'(phi) dphi]
  function pli_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: pli_efold_primitive

    if (alpha.eq.0._kp) stop 'pli_efold_primitive: alpha=0!'

    pli_efold_primitive = -x/alpha

  end function pli_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function pli_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: pli_x_trajectory
   
    
    pli_x_trajectory = xend+alpha*bfold
       
  end function pli_x_trajectory

  
end module plisr
