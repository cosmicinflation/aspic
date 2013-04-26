!slow-roll functions for the inverse monomial potential
!
!V(phi) = M^4 x^(-p)
!
!x = phi/Mp

module imisr
  use infprec, only : kp
  implicit none

  private

  public  imi_norm_potential, imi_epsilon_one, imi_epsilon_two
  public  imi_epsilon_three, imi_x_epsoneunity, imi_xendmin
  public  imi_efold_primitive, imi_x_trajectory
  public  imi_norm_deriv_potential, imi_norm_deriv_second_potential
 
contains
!returns V/M^4
  function imi_norm_potential(x,p)
    implicit none
    real(kp) :: imi_norm_potential
    real(kp), intent(in) :: x,p

    imi_norm_potential = x**(-p)

  end function imi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function imi_norm_deriv_potential(x,p)
    implicit none
    real(kp) :: imi_norm_deriv_potential
    real(kp), intent(in) :: x,p

   imi_norm_deriv_potential = -p*x**(-p-1._kp)

  end function imi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function imi_norm_deriv_second_potential(x,p)
    implicit none
    real(kp) :: imi_norm_deriv_second_potential
    real(kp), intent(in) :: x,p

    imi_norm_deriv_second_potential =  p*(p+1._kp)*x**(-p-2._kp)

  end function imi_norm_deriv_second_potential


!epsilon_one(x)
  function imi_epsilon_one(x,p)    
    implicit none
    real(kp) :: imi_epsilon_one
    real(kp), intent(in) :: x,p
    
    imi_epsilon_one = 0.5_kp*(p/x)**2
    
  end function imi_epsilon_one


!epsilon_two(x)
  function imi_epsilon_two(x,p)    
    implicit none
    real(kp) :: imi_epsilon_two
    real(kp), intent(in) :: x,p
    
    imi_epsilon_two = -2._kp*p/x**2
    
  end function imi_epsilon_two


!epsilon_three(x)
  function imi_epsilon_three(x,p)    
    implicit none
    real(kp) :: imi_epsilon_three
    real(kp), intent(in) :: x,p
    
    imi_epsilon_three = imi_epsilon_two(x,p)
    
  end function imi_epsilon_three


!returns the value of x such that epsilon1=1
  function imi_x_epsoneunity(p)
    implicit none
    real(kp), intent(in) :: p
    real(kp) :: imi_x_epsoneunity
   
    imi_x_epsoneunity = p/sqrt(2._kp)
   
  end function imi_x_epsoneunity


!this is integral[V(phi)/V'(phi) dphi]
  function imi_efold_primitive(x,p)
    implicit none
    real(kp), intent(in) :: x,p
    real(kp) :: imi_efold_primitive

    if (p.eq.0._kp) stop 'imi_efold_primitive: p=0!'

    imi_efold_primitive = -0.5_kp*x**2/p

  end function imi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function imi_x_trajectory(bfold,xend,p)
    implicit none
    real(kp), intent(in) :: bfold, p, xend
    real(kp) :: imi_x_trajectory
           
    imi_x_trajectory = sqrt(2._kp*p*bfold + xend**2)
    
   
  end function imi_x_trajectory

!Returns the minimum value of xend in order to realize the required -bdolstar e-folds.
  function imi_xendmin(efold, p)
    implicit none
    real(kp), intent(in) :: efold, p
    real(kp) :: imi_xendmin

    imi_xendmin=sqrt(0.5_kp*p**2+2._kp*p*efold)

  end function imi_xendmin 

  
end module imisr
