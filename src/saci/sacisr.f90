!slow-roll functions for the SuperConformal alpha-attractor C inflation potential
!
!V(phi) propto [ tanh(x/sqrt(6 alpha)) / ( 1 + tanh(x/sqrt(6 alpha)) )   ]^(2n)
!       = M^4 ( 1 - exp{-sqrt[2/(3 alpha)]x} )^(2n)
!x = phi/Mp

module sacisr
  use infprec, only : kp, tolkp
  use specialinf, only : lambert
  implicit none

  private

  public saci_norm_potential, saci_norm_deriv_potential, saci_norm_deriv_second_potential
  public saci_epsilon_one, saci_epsilon_two, saci_epsilon_three
  public saci_efold_primitive, saci_x_trajectory, saci_x_endinf

  real(kp), parameter :: sqtwothird = sqrt(2._kp/3._kp)
  
contains

!returns V/M**4
  function saci_norm_potential(x,alpha,n)
    implicit none
    real(kp) :: saci_norm_potential,expmy
    real(kp), intent(in) :: x,alpha,n

    expmy = exp(-sqtwothird * x/sqrt(alpha))
    
    saci_norm_potential = (1._kp - expmy)**(2._kp*n)

    
  end function saci_norm_potential


!returns the first derivative of the potential with respect to x=phi/Mp, divided by M**4
  function saci_norm_deriv_potential(x,alpha,n)
    implicit none
    real(kp) :: saci_norm_deriv_potential,expmy
    real(kp), intent(in) :: x,alpha,n

    expmy = exp(-sqtwothird*x/sqrt(alpha))
    
    saci_norm_deriv_potential = 2._kp*n*sqtwothird/sqrt(alpha)*expmy*(1._kp - expmy)**(2._kp*n-1._kp)
    
  end function saci_norm_deriv_potential


!returns the second derivative of the potential with respect to x=phi/Mp, divided by M**4
  function saci_norm_deriv_second_potential(x,alpha,n)
    implicit none
    real(kp) :: saci_norm_deriv_second_potential,expmy,exppy
    real(kp), intent(in) :: x,alpha,n

    expmy = exp(-sqtwothird*x/sqrt(alpha))
    exppy = exp(sqtwothird*x/sqrt(alpha))
    
    saci_norm_deriv_second_potential = 4._kp*n/(3._kp*alpha)*(1._kp - expmy)**(2._kp*n) &
         * (exppy - 2._kp*n) / (exppy - 1._kp)**2
    
  end function saci_norm_deriv_second_potential

!epsilon1(x)
  function saci_epsilon_one(x,alpha,n)
    implicit none
    real(kp) :: saci_epsilon_one,exppy
    real(kp), intent(in) :: x,alpha,n

    exppy = exp(sqtwothird*x/sqrt(alpha))
    
    saci_epsilon_one = (4._kp*n**2)/(3._kp*alpha*(-1._kp + exppy)**2)

  end function saci_epsilon_one


!epsilon2(x)
  function saci_epsilon_two(x,alpha,n)
    implicit none
    real(kp) :: saci_epsilon_two
    real(kp), intent(in) :: x,alpha,n

    saci_epsilon_two = (2._kp*n*1._kp/sinh(x/(sqrt(6._kp*alpha)))**2)/(3._kp*alpha)

  end function saci_epsilon_two

!epsilon3(x)
  function saci_epsilon_three(x,alpha,n)
    implicit none
    real(kp) :: saci_epsilon_three
    real(kp), intent(in) :: x,alpha,n

    saci_epsilon_three = 4._kp*n/(3._kp*alpha)/(exp(sqtwothird*x/sqrt(alpha))-1._kp) &
         /tanh(x/sqrt(6._kp*alpha))

  end function saci_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function saci_x_endinf(alpha,n)
    implicit none
    real(kp) :: saci_x_endinf
    real(kp), intent(in) :: alpha,n

    saci_x_endinf = sqrt(3._kp*alpha/2._kp) * log(1._kp+2._kp*n/sqrt(3._kp*alpha))
    
  end function saci_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function saci_efold_primitive(x,alpha,n)
    implicit none
    real(kp), intent(in) :: x,alpha,n
    real(kp) :: saci_efold_primitive

    saci_efold_primitive = 3._kp*alpha/(4._kp*n)*exp(sqrt(2._kp/(3._kp*alpha))*x) &
                            - sqrt(6._kp*alpha)/(4._kp*n)*x

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
