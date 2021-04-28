!slow-roll functions for the Starobinsky inflation potential
!
!V(phi) = M^4 * [1 - exp(-sqrt(2/3) x)]^2
!
!x = phi/Mp

module sisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  implicit none

  private

  public  si_norm_potential, si_epsilon_one, si_epsilon_two, si_epsilon_three
  public  si_x_endinf, si_efold_primitive, si_x_trajectory
  public  si_norm_deriv_potential, si_norm_deriv_second_potential
 
contains
!returns V/M^4
  function si_norm_potential(x)
    implicit none
    real(kp) :: si_norm_potential
    real(kp), intent(in) :: x

    si_norm_potential = (1._kp-exp(-2._kp*x/(sqrt(6._kp))))**2

  end function si_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function si_norm_deriv_potential(x)
    implicit none
    real(kp) :: si_norm_deriv_potential
    real(kp), intent(in) :: x

   si_norm_deriv_potential = 2._kp*sqrt(2._kp/3._kp)* &
              exp(-2._kp*sqrt(2._kp/3._kp)*x)* (-1._kp+exp(sqrt(2._kp/3._kp)*x))

  end function si_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function si_norm_deriv_second_potential(x)
    implicit none
    real(kp) :: si_norm_deriv_second_potential
    real(kp), intent(in) :: x

    si_norm_deriv_second_potential = -(4._kp/3._kp)*exp(-2._kp*sqrt(2._kp/3._kp)*x)* &
                                      (-2._kp+exp(sqrt(2._kp/3._kp)*x))

  end function si_norm_deriv_second_potential



!epsilon_one(x)
  function si_epsilon_one(x)    
    implicit none
    real(kp) :: si_epsilon_one
    real(kp), intent(in) :: x

  
    si_epsilon_one = 4._kp/(3._kp*(-1._kp+exp(sqrt(2._kp/3._kp)*x))**2)

    
  end function si_epsilon_one


!epsilon_two(x)
  function si_epsilon_two(x)    
    implicit none
    real(kp) :: si_epsilon_two
    real(kp), intent(in) :: x
    
    si_epsilon_two = 2._kp/3._kp/(sinh(x/sqrt(6._kp)))**2
    
  end function si_epsilon_two


!epsilon_three(x)
  function si_epsilon_three(x)    
    implicit none
    real(kp) :: si_epsilon_three
    real(kp), intent(in) :: x
    
    si_epsilon_three = 2._kp/3._kp*(-1+1._kp/tanh(x/sqrt(6._kp)))/tanh(x/sqrt(6._kp))
    
  end function si_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function si_x_endinf()
    implicit none
    real(kp) :: si_x_endinf
    
    si_x_endinf = sqrt(3._kp/2._kp)*log(1._kp+2._kp/sqrt(3._kp))


  end function si_x_endinf



!this is integral[V(phi)/V'(phi) dphi]
  function si_efold_primitive(x)
    implicit none
    real(kp), intent(in) :: x
    real(kp) :: si_efold_primitive

    si_efold_primitive = -0.5_kp*sqrt(1.5_kp)*x+0.75_kp*exp(sqrt(2._kp/3._kp)*x)

  end function si_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function si_x_trajectory(bfold,xend)
    implicit none
    real(kp), intent(in) :: bfold, xend
    real(kp) :: si_x_trajectory
    
    si_x_trajectory = sqrt(1.5_kp)*(4._kp/3._kp*bfold+sqrt(2._kp/3._kp)*xend &
                      -exp(sqrt(2._kp/3._kp)*xend)-lambert(-exp(4._kp/3._kp*bfold &
                      +sqrt(2._kp/3._kp)*xend-exp(sqrt(2._kp/3._kp)*xend)),-1))
       
  end function si_x_trajectory

  



end module sisr
