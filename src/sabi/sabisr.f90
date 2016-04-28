!slow-roll functions for the SuperConformal alpha-attractor B inflation potential
!
!V(phi) =M^4 tanh^(2n) ( x/sqrt(6 alpha) )
!
!x = phi/Mp

module sabisr
  use infprec, only : kp, tolkp
  implicit none

  private

  public sabi_norm_potential, sabi_norm_deriv_potential, sabi_norm_deriv_second_potential
  public sabi_epsilon_one, sabi_epsilon_two, sabi_epsilon_three
  public sabi_efold_primitive, sabi_x_trajectory, sabi_x_endinf

contains

!returns V/M**4
  function sabi_norm_potential(x,alpha,n)
    implicit none
    real(kp) :: sabi_norm_potential
    real(kp), intent(in) :: x,alpha,n

    sabi_norm_potential = tanh(x/sqrt(6._kp*alpha))**(2._kp*n)

  end function sabi_norm_potential


!returns the first derivative of the potential with respect to x=phi/Mp, divided by M**4
  function sabi_norm_deriv_potential(x,alpha,n)
    implicit none
    real(kp) :: sabi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,n

    sabi_norm_deriv_potential = (sqrt(2._kp/3._kp)*n*1._kp/ &
                                cosh(x/(sqrt(6._kp)*sqrt(alpha)))**2*tanh(x/(sqrt(6._kp)* &
                                sqrt(alpha)))**(-1._kp+2._kp*n))/sqrt(alpha)

  end function sabi_norm_deriv_potential


!returns the second derivative of the potential with respect to x=phi/Mp, divided by M**4
  function sabi_norm_deriv_second_potential(x,alpha,n)
    implicit none
    real(kp) :: sabi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,n

sabi_norm_deriv_second_potential = (-4._kp*n*(-2._kp*n+cosh((sqrt(2._kp/3._kp)*x)/sqrt(alpha)))* &
                                1._kp/sinh((sqrt(2._kp/3._kp)*x)/sqrt(alpha))**2*tanh(x/ &
                                (sqrt(6._kp)*sqrt(alpha)))**(2._kp*n))/(3._kp*alpha)

  end function sabi_norm_deriv_second_potential

!epsilon1(x)
  function sabi_epsilon_one(x,alpha,n)
    implicit none
    real(kp) :: sabi_epsilon_one
    real(kp), intent(in) :: x,alpha,n

    sabi_epsilon_one = (4._kp*n**2*1._kp/sinh((sqrt(2._kp/3._kp)*x)/sqrt(alpha))**2)/(3._kp*alpha)

  end function sabi_epsilon_one


!epsilon2(x)
  function sabi_epsilon_two(x,alpha,n)
    implicit none
    real(kp) :: sabi_epsilon_two
    real(kp), intent(in) :: x,alpha,n

    sabi_epsilon_two = (2._kp*n*(2._kp+1._kp/sinh(x/(sqrt(6._kp)*sqrt(alpha)))**2)* &
                        1._kp/cosh(x/(sqrt(6._kp)*sqrt(alpha)))**2)/(3._kp*alpha)

  end function sabi_epsilon_two

!epsilon3(x)
  function sabi_epsilon_three(x,alpha,n)
    implicit none
    real(kp) :: sabi_epsilon_three
    real(kp), intent(in) :: x,alpha,n

    sabi_epsilon_three = (2._kp*n*(3._kp+cosh((2._kp*sqrt(2._kp/3._kp)*x)/sqrt(alpha)))* &
                        1._kp/sinh((sqrt(2._kp/3._kp)*x)/sqrt(alpha))**2*1._kp/cosh((sqrt(2._kp/ &
                        3._kp)*x)/sqrt(alpha)))/(3._kp*alpha)

  end function sabi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function sabi_x_endinf(alpha,n)
    implicit none
    real(kp) :: sabi_x_endinf
    real(kp), intent(in) :: alpha,n

    sabi_x_endinf = sqrt(6._kp*alpha)/2._kp*asinh(2._kp*n/sqrt(3._kp*alpha))

  end function sabi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function sabi_efold_primitive(x,alpha,n)
    implicit none
    real(kp), intent(in) :: x,alpha,n
    real(kp) :: sabi_efold_primitive

    sabi_efold_primitive = 3._kp*alpha/(4._kp*n)*cosh(sqrt(2._kp/(3._kp*alpha))*x)

  end function sabi_efold_primitive

!returns x at bfold=-efolds before the end of inflation
  function sabi_x_trajectory(bfold,xend,alpha,n)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,n
    real(kp) :: sabi_x_trajectory

    sabi_x_trajectory = sqrt(3._kp*alpha/2._kp)*acosh(sqrt(1._kp+4._kp*n**2/alpha)- &
                        4._kp*bfold*n/(3._kp*alpha))

  end function sabi_x_trajectory


end module sabisr
