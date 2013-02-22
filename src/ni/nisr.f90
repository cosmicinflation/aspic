!slow-roll functions for the natural inflation potential
!
!V(phi) = M^4 [ 1 + cos(x/f) ]
!
!x = phi/Mp
!f = f/Mp

module nisr
  use infprec, only : kp
  implicit none

  private

  public  ni_norm_potential, ni_epsilon_one, ni_epsilon_two, ni_epsilon_three
  public  ni_x_endinf, ni_efold_primitive, ni_x_trajectory
  public  ni_norm_deriv_potential, ni_norm_deriv_second_potential
 
contains
!returns V/M^4
  function ni_norm_potential(x,f)
    implicit none
    real(kp) :: ni_norm_potential
    real(kp), intent(in) :: x,f

    ni_norm_potential = 1._kp+cos(x/f)

  end function ni_norm_potential

!returns the first derivative of the potential with respect to x, divided by M^4
  function ni_norm_deriv_potential(x,f)
    implicit none
    real(kp) :: ni_norm_deriv_potential
    real(kp), intent(in) :: x,f

    ni_norm_deriv_potential = -sin(x/f)/f

  end function ni_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function ni_norm_deriv_second_potential(x,f)
    implicit none
    real(kp) :: ni_norm_deriv_second_potential
    real(kp), intent(in) :: x,f

    ni_norm_deriv_second_potential = -cos(x/f)/(f**2)

  end function ni_norm_deriv_second_potential


!epsilon_one(x)
  function ni_epsilon_one(x,f)    
    implicit none
    real(kp) :: ni_epsilon_one
    real(kp), intent(in) :: x,f
    
    ni_epsilon_one = 1._kp/(2._kp*f**2) &
         *(sin(x/f))**2/(1._kp+cos(x/f))**2
    
  end function ni_epsilon_one


!epsilon_two(x)
  function ni_epsilon_two(x,f)
    implicit none
    real(kp) :: ni_epsilon_two
    real(kp), intent(in) :: x,f
    
    ni_epsilon_two = 2._kp/(f**2) &
         /(1._kp+cos(x/f))
    
  end function ni_epsilon_two


!epsilon_three(x)
  function ni_epsilon_three(x,f)
    implicit none
    real(kp) :: ni_epsilon_three
    real(kp), intent(in) :: x,f
    
    ni_epsilon_three = 2._kp*ni_epsilon_one(x,f) 
    
  end function ni_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function ni_x_endinf(f)
    implicit none
    real(kp), intent(in) :: f
    real(kp) :: ni_x_endinf
   
    ni_x_endinf = f*acos((1._kp-2._kp*f**2)/(1._kp+2._kp*f**2))
   
  end function ni_x_endinf



!this is integral[V(phi)/V'(phi) dphi]
  function ni_efold_primitive(x,f)
    implicit none
    real(kp), intent(in) :: x,f
    real(kp) :: ni_efold_primitive
    
    ni_efold_primitive = -f**2 * log(abs( &
         1._kp-cos(x/f) ))
         
   
  end function ni_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ni_x_trajectory(bfold,xend,f)
    implicit none
    real(kp), intent(in) :: bfold,f,xend
    real(kp) :: ni_x_trajectory
    
    ni_x_trajectory = f*acos(1._kp-(1._kp-cos(xend/f))*exp(bfold/(f**2)))
       
  end function ni_x_trajectory

  
  
end module nisr
