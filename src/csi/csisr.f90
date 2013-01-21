!slow-roll functions for the constant spectrum inflation potential
!
!V(phi) = M**4 / (1 - alpha phi/Mp)**2
!
!x = phi/Mp

module csisr
  use infprec, only : kp, tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public csi_norm_potential, csi_norm_deriv_potential, csi_norm_deriv_second_potential
  public csi_epsilon_one, csi_epsilon_two,csi_epsilon_three
  public csi_efold_primitive, csi_x_trajectory
  public csi_xendmax, csi_x_epsoneunity

 
contains
!returns V/M**4
  function csi_norm_potential(x,alpha)
    implicit none
    real(kp) :: csi_norm_potential
    real(kp), intent(in) :: x,alpha

    csi_norm_potential = 1._kp/(1._kp-alpha*x)**2
  end function csi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M**4
  function csi_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: csi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

   csi_norm_deriv_potential = (2._kp*alpha)/(1._kp-alpha*x)**3

  end function csi_norm_deriv_potential


!returns the second derivative of the potential with respect to x, divided by M**4
  function csi_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: csi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    csi_norm_deriv_second_potential = (6._kp*alpha**2)/(1._kp-alpha*x)**4

  end function csi_norm_deriv_second_potential

!epsilon1(x)
  function csi_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: csi_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    csi_epsilon_one = (2._kp*alpha**2)/(-1._kp+alpha*x)**2
    
  end function csi_epsilon_one


!epsilon2(x)
  function csi_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: csi_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    csi_epsilon_two = -((4._kp*alpha**2)/(-1._kp+alpha*x)**2)
    
  end function csi_epsilon_two

!epsilon3(x)
  function csi_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: csi_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    csi_epsilon_three = -((4._kp*alpha**2)/(-1._kp+alpha*x)**2)
    
  end function csi_epsilon_three


!this is integral[V(phi)/V'(phi) dphi]
  function csi_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: csi_efold_primitive
    
     if (alpha.eq.0._kp) stop 'csi_efold_primitive: alpha=0 is singular'

    csi_efold_primitive = x/(2._kp*alpha)-x**2/4._kp

  end function csi_efold_primitive
 


!returns x at bfold=-efolds before the end of inflation
  function csi_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha
    real(kp) :: csi_x_trajectory
    
    csi_x_trajectory = (1._kp-sqrt(1._kp-2._kp*alpha*xend+ &
                       alpha**2*xend**2+4._kp*alpha**2*bfold))/alpha !in the x<1/alpha branch of the potential

!    csi_x_trajectory = (1._kp+sqrt(1._kp-2._kp*alpha*xend+ &
!                       alpha**2*xend**2+4._kp*alpha**2*bfold))/alpha !in the x>1/alpha branch of the potential
    
  end function csi_x_trajectory

! Returns the value of x such that epsilon_one=1
  function csi_x_epsoneunity(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: csi_x_epsoneunity

    csi_x_epsoneunity=1/alpha-sqrt(2._kp) !in the x<1/alpha branch of the potential

!    csi_x_epsOne_equals_one=1/alpha+sqrt(2._kp) !in the x>1/alpha branch of the potential

  end function csi_x_epsoneunity

! Returns the value of x such that epsilon_two=epsilon_three=-1
  function csi_x_epsTwo_equals_MinusOne(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: csi_x_epsTwo_equals_MinusOne

    csi_x_epsTwo_equals_MinusOne=1/alpha-2._kp !in the x<1/alpha branch of the potential

!    csi_x_epsTwo_equals_MinusOne=1/alpha+2._kp !in the x>1/alpha branch of the potential

  end function csi_x_epsTwo_equals_MinusOne

  
  function csi_xendmax(efold,alpha) !Returns the maximum value of xend in order to realize the required -bdolstar e-folds.
    implicit none
    real(kp), intent(in) :: alpha, efold
    real(kp) :: csi_xendmax
    
    csi_xendmax = csi_x_trajectory(efold,csi_x_epsoneunity(alpha),alpha)

  end function csi_xendmax



end module csisr
