!slow-roll functions for the beta exponential inflation potential
!
!V(phi) = M^4 [ 1 - beta lambda x ] ^ ( 1/beta )
!
!x = phi/Mp

module beisr
  use infprec, only : kp, tolkp,transfert
  implicit none

  private

  public bei_norm_potential, bei_norm_deriv_potential, bei_norm_deriv_second_potential
  public bei_epsilon_one, bei_epsilon_two,bei_epsilon_three
  public bei_efold_primitive, bei_x_trajectory,bei_x_endinf

 
contains
!returns V/M**4
  function bei_norm_potential(x,lambda,beta)
    implicit none
    real(kp) :: bei_norm_potential
    real(kp), intent(in) :: x,lambda,beta

    bei_norm_potential = (1._kp-beta*lambda*x)**(1._kp/beta)
  end function bei_norm_potential


!returns the first derivative of the potential with respect to x=phi/Mp, divided by M**4
  function bei_norm_deriv_potential(x,lambda,beta)
    implicit none
    real(kp) :: bei_norm_deriv_potential
    real(kp), intent(in) :: x,lambda,beta

   bei_norm_deriv_potential = -lambda*(1._kp-beta*lambda*x)**(1._kp/beta-1._kp)

  end function bei_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function bei_norm_deriv_second_potential(x,lambda,beta)
    implicit none
    real(kp) :: bei_norm_deriv_second_potential
    real(kp), intent(in) :: x,lambda,beta

    bei_norm_deriv_second_potential = (1._kp-beta)*lambda**2* &
                                      (1._kp-beta*lambda*x)**(1._kp/beta-2._kp)

  end function bei_norm_deriv_second_potential

!epsilon1(x)
  function bei_epsilon_one(x,lambda,beta)    
    implicit none
    real(kp) :: bei_epsilon_one
    real(kp), intent(in) :: x,lambda,beta
    
    bei_epsilon_one = lambda**2/(2._kp*(1._kp-beta*lambda*x)**2)
    
  end function bei_epsilon_one


!epsilon2(x)
  function bei_epsilon_two(x,lambda,beta)    
    implicit none
    real(kp) :: bei_epsilon_two
    real(kp), intent(in) :: x,lambda,beta
    
    bei_epsilon_two = 4._kp*beta*bei_epsilon_one(x,lambda,beta) 
    
  end function bei_epsilon_two

!epsilon3(x)
  function bei_epsilon_three(x,lambda,beta)    
    implicit none
    real(kp) :: bei_epsilon_three
    real(kp), intent(in) :: x,lambda,beta
    
    bei_epsilon_three = bei_epsilon_two(x,lambda,beta)    
    
  end function bei_epsilon_three

!returns x at the end of inflation defined as epsilon1=1
  function bei_x_endinf(lambda,beta)
    implicit none
    real(kp) :: bei_x_endinf
    real(kp), intent(in) :: lambda,beta
   
     if (beta.eq.0._kp) stop 'bei_x_endinf: beta=0 is singular'

    bei_x_endinf = (1._kp/lambda-1._kp/sqrt(2._kp))/beta

  end function bei_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function bei_efold_primitive(x,lambda,beta)
    implicit none
    real(kp), intent(in) :: x,lambda,beta
    real(kp) :: bei_efold_primitive
    
     if (lambda.eq.0._kp) stop 'bei_efold_primitive: lambda=0 is singular'

    bei_efold_primitive = -x/lambda+0.5_kp*beta*x**2

  end function bei_efold_primitive
 


!returns x at bfold=-efolds before the end of inflation
  function bei_x_trajectory(bfold,xend,lambda,beta)
    implicit none
    real(kp), intent(in) :: bfold,lambda,beta,xend
    real(kp) :: bei_x_trajectory
    
    
    bei_x_trajectory = 1._kp/lambda/beta - sqrt((xend - 1._kp/lambda/beta)**2 - 2._kp*bfold/beta)
    
  end function bei_x_trajectory



end module beisr
