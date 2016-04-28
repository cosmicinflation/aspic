!slow-roll functions for the SuperConformal alpha-attractor C inflation potential
!
!V(phi) =  M^4 [ tanh(x/sqrt(6 alpha)) / ( 1 + tanh(x/sqrt(6 alpha)) )   ]^(2n)
!
!x = phi/Mp

module sacisr
  use infprec, only : kp, tolkp
  use specialinf, only : lambert
  implicit none

  private

  public saci_norm_potential, saci_norm_deriv_potential, saci_norm_deriv_second_potential
  public saci_epsilon_one, saci_epsilon_two, saci_epsilon_three
  public saci_efold_primitive, saci_x_trajectory, saci_x_endinf

contains

!returns V/M**4
  function saci_norm_potential(x,alpha,n)
    implicit none
    real(kp) :: saci_norm_potential, y
    real(kp), intent(in) :: x,alpha,n

    y = tanh(x/sqrt(6._kp*alpha))
    saci_norm_potential = (y/(y+1._kp))**(2._kp*n)

  end function saci_norm_potential


!returns the first derivative of the potential with respect to x=phi/Mp, divided by M**4
  function saci_norm_deriv_potential(x,alpha,n)
    implicit none
    real(kp) :: saci_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,n

    saci_norm_deriv_potential = (sqrt(2._kp/3._kp)*n*1._kp/cosh(x/(sqrt(6._kp)*sqrt(alpha)))**2* &
                                (tanh(x/(sqrt(6._kp)*sqrt(alpha)))/(1._kp+tanh(x/(sqrt(6._kp)* &
                                sqrt(alpha)))))**(-1._kp+2._kp*n))/(sqrt(alpha)*(1._kp+ &
                                tanh(x/(sqrt(6._kp)*sqrt(alpha))))**2)

  end function saci_norm_deriv_potential


!returns the second derivative of the potential with respect to x=phi/Mp, divided by M**4
  function saci_norm_deriv_second_potential(x,alpha,n)
    implicit none
    real(kp) :: saci_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,n

    saci_norm_deriv_second_potential = (n*1._kp/sinh(x/(sqrt(6._kp)*sqrt(alpha)))**2*(tanh(x/ &
                                        (sqrt(6._kp)*sqrt(alpha)))/(1._kp+tanh(x/(sqrt(6._kp)* &
                                        sqrt(alpha)))))**(2*n)*(2._kp*n*1._kp/cosh(x/ &
                                        (sqrt(6._kp)*sqrt(alpha)))**2-(1._kp+tanh(x/(sqrt(6._kp)* &
                                        sqrt(alpha))))**2))/(3._kp*alpha*(1._kp+tanh(x/(sqrt(6._kp)* &
                                        sqrt(alpha))))**2)

  end function saci_norm_deriv_second_potential

!epsilon1(x)
  function saci_epsilon_one(x,alpha,n)
    implicit none
    real(kp) :: saci_epsilon_one
    real(kp), intent(in) :: x,alpha,n

    saci_epsilon_one = (4._kp*n**2)/(3._kp*alpha*(-1._kp+exp((sqrt(2._kp/3._kp)*x)/sqrt(alpha)))**2)

  end function saci_epsilon_one


!epsilon2(x)
  function saci_epsilon_two(x,alpha,n)
    implicit none
    real(kp) :: saci_epsilon_two
    real(kp), intent(in) :: x,alpha,n

    saci_epsilon_two = (2._kp*n*1._kp/sinh(x/(sqrt(6._kp*alpha)))**2)/(3.*alpha)

  end function saci_epsilon_two

!epsilon3(x)
  function saci_epsilon_three(x,alpha,n)
    implicit none
    real(kp) :: saci_epsilon_three
    real(kp), intent(in) :: x,alpha,n

    saci_epsilon_three = (2*n*(-1._kp+1._kp/tanh(x/(sqrt(6._kp)*sqrt(alpha))))* &
                            1._kp/tanh(x/(sqrt(6._kp)*sqrt(alpha))))/(3._kp*alpha)

  end function saci_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function saci_x_endinf(alpha,n)
    implicit none
    real(kp) :: saci_x_endinf
    real(kp), intent(in) :: alpha,n

    saci_x_endinf = sqrt(6._kp*alpha)/2._kp*log(1._kp+2._kp*n/sqrt(3._kp*alpha))
  end function saci_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function saci_efold_primitive(x,alpha,n)
    implicit none
    real(kp), intent(in) :: x,alpha,n
    real(kp) :: saci_efold_primitive

    saci_efold_primitive = 3._kp*alpha/(4._kp*n)*exp(sqrt(2._kp/(3._kp*alpha))*x)- &
                            sqrt(6._kp*alpha)/(4._kp*n)*x

  end function saci_efold_primitive

!returns x at bfold=-efolds before the end of inflation
  function saci_x_trajectory(bfold,xend,alpha,n)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,n
    real(kp) :: saci_x_trajectory

    saci_x_trajectory = -(sqrt(2._kp)*n)+(2._kp*sqrt(2._kp/3._kp)*bfold*n)/sqrt(alpha)+ &
                        sqrt(1.5)*sqrt(alpha)*(-1._kp+log(1._kp+(2._kp*n)/(sqrt(3._kp)* &
                        sqrt(alpha)))-lambert((exp(-1._kp-(2._kp*(sqrt(3._kp)*sqrt(alpha)- &
                        2._kp*bfold)*n)/(3._kp*alpha))*(-3._kp-(2._kp*sqrt(3._kp)*n)/ &
                        sqrt(alpha)))/3._kp,-1))

  end function saci_x_trajectory


end module sacisr
