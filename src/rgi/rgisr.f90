!slow-roll functions for the radion gauge potential
!
!V(phi) = M^4 x^2 / [ alpha + x^2 ]
!
!x = phi/Mp

module rgisr
  use infprec, only : kp
  implicit none

  private

  public rgi_norm_potential, rgi_epsilon_one, rgi_epsilon_two
  public rgi_epsilon_three
  public rgi_x_endinf, rgi_efold_primitive, rgi_x_trajectory
  public rgi_norm_deriv_potential, rgi_norm_deriv_second_potential
  public rgi_x_epsonenumacc

contains
!returns V/M^4
  function rgi_norm_potential(x,alpha)
    implicit none
    real(kp) :: rgi_norm_potential
    real(kp), intent(in) :: x,alpha

    rgi_norm_potential = x**2/(alpha+x**2)

  end function rgi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function rgi_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: rgi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

   rgi_norm_deriv_potential = 2._kp*alpha*x/((alpha+x**2)**2)

  end function rgi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function rgi_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: rgi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    rgi_norm_deriv_second_potential =  2._kp*alpha*(alpha-3._kp*x**2) &
         /((alpha+x**2)**3)

  end function rgi_norm_deriv_second_potential


!epsilon_one(x)
  function rgi_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: rgi_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    rgi_epsilon_one = 2._kp*alpha**2/(x**2*(alpha+x**2)**2)
    
  end function rgi_epsilon_one


!epsilon_two(x)
  function rgi_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: rgi_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    rgi_epsilon_two = 4._kp*alpha*(alpha+3._kp*x**2)/(x**2*(alpha+x**2)**2)
    
  end function rgi_epsilon_two


!epsilon_three(x)
  function rgi_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: rgi_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    rgi_epsilon_three = 4._kp*alpha*(alpha**2+3._kp*alpha*x**2+6._kp*x**4)/ &
         (x**2*(alpha+x**2)**2*(alpha+3._kp*x**2))
    
  end function rgi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function rgi_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: rgi_x_endinf
   
    rgi_x_endinf = (6._kp**(1._kp/3._kp)*alpha-(alpha*(-9._kp+ &
         sqrt(81._kp+6._kp*alpha)))**(2._kp/3._kp))/ &
         (2._kp**(1._kp/6._kp)*3**(2._kp/3_kp)*(alpha* &
         (-9._kp+sqrt(81._kp+6_kp*alpha)))**(1._kp/3._kp))
    
  end function rgi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function rgi_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: rgi_efold_primitive

    if (alpha.eq.0._kp) stop 'rgi_efold_primitive: alpha=0!'

    rgi_efold_primitive = x**2/4._kp+x**4/(8._kp*alpha)

  end function rgi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function rgi_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: rgi_x_trajectory
           
    rgi_x_trajectory = sqrt(-alpha+sqrt(alpha**2+xend**4+ &
                       2._kp*alpha*(-4._kp*bfold+xend**2)))
    
   
  end function rgi_x_trajectory


!return x such that epsilon_1 = 1e-16 for checking numerical accuracy  
  function rgi_x_epsonenumacc(alpha)
    implicit none
    real(kp) :: rgi_x_epsonenumacc
    real(kp), intent(in) :: alpha
    real(kp) :: epsMini

    epsMini = 10._kp*epsilon(1._kp)

    rgi_x_epsonenumacc = (2._kp*alpha*alpha/epsMini)**(1._kp/6._kp)
    
  end function rgi_x_epsonenumacc
  
end module rgisr
