!slow-roll functions for the SuperConformal alpha Attractor B inflation potential
!
!V(phi) =  M^4 * ( 1 - exp(-sqrt(2/(3 alpha)) x ) )^2
!
!x = phi/Mp
!
!alpha is a dimensionless parameter that is 1 for Starobinsky Inflation

module saaisr
  use infprec, only : kp, tolkp, transfert
  use specialinf, only : lambert
  implicit none

  private

  public  saai_norm_potential, saai_epsilon_one, saai_epsilon_two, saai_epsilon_three
  public  saai_x_endinf, saai_efold_primitive, saai_x_trajectory
  public  saai_norm_deriv_potential, saai_norm_deriv_second_potential


contains


!returns V/M**4
  function saai_norm_potential(x, alpha)
    implicit none
    real(kp) :: saai_norm_potential
    real(kp), intent(in) :: x, alpha

    saai_norm_potential = (1._kp-exp(-sqrt(2._kp/(3._kp*alpha))*x))**2

  end function saai_norm_potential


!returns the first derivative of the potential with respect to x, divided by M**4
  function saai_norm_deriv_potential(x, alpha)
    implicit none
    real(kp) :: saai_norm_deriv_potential
    real(kp), intent(in) :: x, alpha

    saai_norm_deriv_potential = (2._kp*sqrt(2._kp/3._kp)*(-1._kp+exp((sqrt(2._kp/3._kp)*x)/ &
                                sqrt(alpha))))/(sqrt(alpha)*exp((2._kp*sqrt(2._kp/3._kp)*x)/ &
                                sqrt(alpha)))

  end function saai_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function saai_norm_deriv_second_potential(x, alpha)
    implicit none
    real(kp) :: saai_norm_deriv_second_potential
    real(kp), intent(in) :: x, alpha

    saai_norm_deriv_second_potential = (-4._kp*(-2._kp+exp((sqrt(2._kp/3._kp)*x)/ &
                sqrt(alpha))))/(3._kp*alpha*exp((2._kp*sqrt(2._kp/3._kp)*x)/sqrt(alpha)))

  end function saai_norm_deriv_second_potential



!epsilon_one(x)
  function saai_epsilon_one(x, alpha)    
    implicit none
    real(kp) :: saai_epsilon_one
    real(kp), intent(in) :: x, alpha
    
    saai_epsilon_one = 4._kp/(3._kp*alpha*(-1._kp+exp(sqrt(2._kp/3._kp)*sqrt(1/alpha)*x))**2)
    
  end function saai_epsilon_one


!epsilon_two(x)
  function saai_epsilon_two(x, alpha)    
    implicit none
    real(kp) :: saai_epsilon_two
    real(kp), intent(in) :: x, alpha
    
    saai_epsilon_two = (2._kp/sinh(x/(sqrt(6._kp)*sqrt(alpha)))**2)/(3._kp*alpha)
    
  end function saai_epsilon_two


!epsilon_three(x)
  function saai_epsilon_three(x, alpha)    
    implicit none
    real(kp) :: saai_epsilon_three
    real(kp), intent(in) :: x, alpha
    
    saai_epsilon_three = (2._kp*(-1._kp+1._kp/tanh((sqrt(1._kp/alpha)*x)/sqrt(6._kp)))* &
                        1._kp/tanh((sqrt(1._kp/alpha)*x)/sqrt(6._kp)))/(3._kp*alpha)
    
  end function saai_epsilon_three

!returns x at the end of inflation defined as epsilon1=1
  function saai_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: saai_x_endinf


    saai_x_endinf = sqrt(3._kp*alpha/2._kp)*log(1._kp+2._kp/sqrt(3._kp*alpha))
   
  end function saai_x_endinf

!this is integral(V(phi)/V'(phi) dphi)
  function saai_efold_primitive(x, alpha)
    implicit none
    real(kp), intent(in) :: x, alpha
    real(kp) :: saai_efold_primitive

    saai_efold_primitive = 0.75_kp*alpha*exp(sqrt(2._kp/(3._kp*alpha))*x)- &
                            sqrt(6._kp*alpha)/4._kp*x

  end function saai_efold_primitive


!returns x at bfold=-efolds before the end of inflation
  function saai_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha
    real(kp) :: saai_x_trajectory, xx

    xx = 4._kp/(3._kp*alpha)*bfold+log(1._kp+2._kp/sqrt(3._kp*alpha))-(1._kp+2._kp/sqrt(3._kp*alpha))

    saai_x_trajectory = sqrt(3._kp*alpha/2._kp)*(xx-lambert(-exp(xx),-1))

  end function saai_x_trajectory

  
end module saaisr
