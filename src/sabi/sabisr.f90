!slow-roll functions for the SuperConformal alpha-attractor C inflation potential
!
!V(phi) = Mo^4 [tanh(x/sqrt(6 alpha)) / ( 1 + tanh(x/sqrt(6 alpha)) )   ]^(2n)
!
!       = M^4 ( 1 - exp{-sqrt[2/(3 alpha)]x} )^(2n)
!
!x = phi/Mp

module sabisr
  use infprec, only : kp, tolkp
  use specialinf, only : lambert
  implicit none

  private

  public sabi_norm_potential, sabi_norm_deriv_potential, sabi_norm_deriv_second_potential
  public sabi_epsilon_one, sabi_epsilon_two, sabi_epsilon_three
  public sabi_efold_primitive, sabi_x_trajectory, sabi_x_endinf

  real(kp), parameter :: sqtwothird = sqrt(2._kp/3._kp)
  
contains

!returns V/M**4
  function sabi_norm_potential(x,alpha,n)
    implicit none
    real(kp) :: sabi_norm_potential,expmy
    real(kp), intent(in) :: x,alpha,n

    expmy = exp(-sqtwothird * x/sqrt(alpha))
    
    sabi_norm_potential = (1._kp - expmy)**(2._kp*n)

    
  end function sabi_norm_potential


!returns the first derivative of the potential with respect to x=phi/Mp, divided by M**4
  function sabi_norm_deriv_potential(x,alpha,n)
    implicit none
    real(kp) :: sabi_norm_deriv_potential,expmy
    real(kp), intent(in) :: x,alpha,n

    expmy = exp(-sqtwothird*x/sqrt(alpha))
    
    sabi_norm_deriv_potential = 2._kp*n*sqtwothird/sqrt(alpha)*expmy*(1._kp - expmy)**(2._kp*n-1._kp)
    
  end function sabi_norm_deriv_potential


!returns the second derivative of the potential with respect to x=phi/Mp, divided by M**4
  function sabi_norm_deriv_second_potential(x,alpha,n)
    implicit none
    real(kp) :: sabi_norm_deriv_second_potential,expmy,exppy
    real(kp), intent(in) :: x,alpha,n

    expmy = exp(-sqtwothird*x/sqrt(alpha))
    exppy = exp(sqtwothird*x/sqrt(alpha))
    
    sabi_norm_deriv_second_potential = -4._kp*n/(3._kp*alpha)*(1._kp - expmy)**(2._kp*n) &
         * (exppy - 2._kp*n) / (exppy - 1._kp)**2
    
  end function sabi_norm_deriv_second_potential

!epsilon1(x)
  function sabi_epsilon_one(x,alpha,n)
    implicit none
    real(kp) :: sabi_epsilon_one,exppy
    real(kp), intent(in) :: x,alpha,n

    exppy = exp(sqtwothird*x/sqrt(alpha))
    
    sabi_epsilon_one = (4._kp*n**2)/(3._kp*alpha*(-1._kp + exppy)**2)

  end function sabi_epsilon_one


!epsilon2(x)
  function sabi_epsilon_two(x,alpha,n)
    implicit none
    real(kp) :: sabi_epsilon_two
    real(kp), intent(in) :: x,alpha,n

    sabi_epsilon_two = (2._kp*n*1._kp/sinh(x/(sqrt(6._kp*alpha)))**2)/(3._kp*alpha)

  end function sabi_epsilon_two

!epsilon3(x)
  function sabi_epsilon_three(x,alpha,n)
    implicit none
    real(kp) :: sabi_epsilon_three
    real(kp), intent(in) :: x,alpha,n

    sabi_epsilon_three = 4._kp*n/(3._kp*alpha)/(exp(sqtwothird*x/sqrt(alpha))-1._kp) &
         /tanh(x/sqrt(6._kp*alpha))

  end function sabi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function sabi_x_endinf(alpha,n)
    implicit none
    real(kp) :: sabi_x_endinf
    real(kp), intent(in) :: alpha,n

    sabi_x_endinf = sqrt(3._kp*alpha/2._kp) * log(1._kp+2._kp*n/sqrt(3._kp*alpha))
    
  end function sabi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function sabi_efold_primitive(x,alpha,n)
    implicit none
    real(kp), intent(in) :: x,alpha,n
    real(kp) :: sabi_efold_primitive

    sabi_efold_primitive = 3._kp*alpha/(4._kp*n)*exp(sqrt(2._kp/(3._kp*alpha))*x) &
         - sqrt(6._kp*alpha)/(4._kp*n)*x

  end function sabi_efold_primitive

!returns x at bfold=-efolds before the end of inflation
  function sabi_x_trajectory(bfold,xend,alpha,n)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,n
    real(kp) :: sabi_x_trajectory

    sabi_x_trajectory = -(sqrt(2._kp)*n)+(2._kp*sqrt(2._kp/3._kp)*bfold*n)/sqrt(alpha)+ &
         sqrt(1.5)*sqrt(alpha)*(-1._kp+log(1._kp+(2._kp*n)/(sqrt(3._kp)* &
         sqrt(alpha)))-lambert((exp(-1._kp-(2._kp*(sqrt(3._kp)*sqrt(alpha)- &
         2._kp*bfold)*n)/(3._kp*alpha))*(-3._kp-(2._kp*sqrt(3._kp)*n)/ &
         sqrt(alpha)))/3._kp,-1))

  end function sabi_x_trajectory


end module sabisr
