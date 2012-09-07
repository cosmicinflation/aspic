!slow-roll functions for the exponential SUSY inflation potential
!
!V(phi) = M^4 * [1- exp(-q * x) ]
!
!x = phi/Mp

module esisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  implicit none

  private

  public  esi_norm_potential, esi_epsilon_one, esi_epsilon_two, esi_epsilon_three
  public  esi_x_endinf, esi_efold_primitive, esi_x_trajectory
  public  esi_norm_deriv_potential, esi_norm_deriv_second_potential
 
contains
!returns V/M^4
  function esi_norm_potential(x,q)
    implicit none
    real(kp) :: esi_norm_potential
    real(kp), intent(in) :: x,q

    esi_norm_potential = 1._kp-exp(-q*x)

  end function esi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function esi_norm_deriv_potential(x,q)
    implicit none
    real(kp) :: esi_norm_deriv_potential
    real(kp), intent(in) :: x,q

   esi_norm_deriv_potential = q*exp(-q*x)

  end function esi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function esi_norm_deriv_second_potential(x,q)
    implicit none
    real(kp) :: esi_norm_deriv_second_potential
    real(kp), intent(in) :: x,q

    esi_norm_deriv_second_potential =  -q**2*exp(-q*x)

  end function esi_norm_deriv_second_potential



!epsilon_one(x)
  function esi_epsilon_one(x,q)    
    implicit none
    real(kp) :: esi_epsilon_one
    real(kp), intent(in) :: x,q
    
    esi_epsilon_one = q**2/2._kp*exp(-2._kp*q*x) &
         /(1._kp-exp(-q*x))**2

    
  end function esi_epsilon_one


!epsilon_two(x)
  function esi_epsilon_two(x,q)    
    implicit none
    real(kp) :: esi_epsilon_two
    real(kp), intent(in) :: x,q
    
    esi_epsilon_two = 2._kp*q**2*exp(-q*x) &
         /(1._kp-exp(-q*x))**2
    
  end function esi_epsilon_two


!epsilon_three(x)
  function esi_epsilon_three(x,q)    
    implicit none
    real(kp) :: esi_epsilon_three
    real(kp), intent(in) :: x,q
    
    esi_epsilon_three =q**2*exp(-q*x) &
         *(1._kp+exp(-q*x)) &
         /(1._kp-exp(-q*x))**2
    
  end function esi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function esi_x_endinf(q)
    implicit none
    real(kp), intent(in) :: q
    real(kp) :: esi_x_endinf
   
    esi_x_endinf = 1._kp/q*log(1._kp+q/sqrt(2._kp))
   
  end function esi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function esi_efold_primitive(x,q)
    implicit none
    real(kp), intent(in) :: x,q
    real(kp) :: esi_efold_primitive

    if (q.eq.0._kp) stop 'esi_efold_primitive: q=0!'

    esi_efold_primitive = 1._kp/q**2 * (exp(q*x)-q*x)

  end function esi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function esi_x_trajectory(bfold,xend,q)
    implicit none
    real(kp), intent(in) :: bfold, q, xend
    real(kp) :: esi_x_trajectory
       
    esi_x_trajectory = q*bfold+xend-exp(q*xend)/q &
         -1._kp/q*lambert(-exp(q**2*bfold+q*xend-exp(q*xend)),-1)
       
  end function esi_x_trajectory



  
end module esisr
