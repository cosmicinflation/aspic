!slow-roll functions for the horizon flow inflation at first order potential
!
!V(phi) = M^4 (1 + A1 phi)^2 * [1 - 2/3 Mp^2 (A1/(1+A1 phi))^2 ]
!
!x = phi/Mp

module hf1isr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : ei
  use inftools, only : zbrent
  implicit none

  private

  public  hf1i_norm_potential, hf1i_epsilon_one, hf1i_epsilon_two, hf1i_epsilon_three
  public  hf1i_x_endinf, hf1i_efold_primitive, hf1i_x_trajectory
  public  hf1i_norm_deriv_potential, hf1i_norm_deriv_second_potential
 
contains
!returns V/M^4
  function hf1i_norm_potential(x,A1)
    implicit none
    real(kp) :: hf1i_norm_potential
    real(kp), intent(in) :: x,A1

    hf1i_norm_potential = (1._kp+A1*x)**2 &
         *(1._kp-2._kp/3._kp*(A1/(1._kp+A1*x))**2)

  end function hf1i_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function hf1i_norm_deriv_potential(x,A1)
    implicit none
    real(kp) :: hf1i_norm_deriv_potential
    real(kp), intent(in) :: x,A1

   hf1i_norm_deriv_potential = 2._kp*A1*(1._kp+A1*x)

  end function hf1i_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function hf1i_norm_deriv_second_potential(x,A1)
    implicit none
    real(kp) :: hf1i_norm_deriv_second_potential
    real(kp), intent(in) :: x,A1

    hf1i_norm_deriv_second_potential =  2._kp*A1**2

  end function hf1i_norm_deriv_second_potential




!epsilon_one(x)
  function hf1i_epsilon_one(x,A1)    
    implicit none
    real(kp) :: hf1i_epsilon_one
    real(kp), intent(in) :: x,A1
    
    hf1i_epsilon_one = 2._kp*(A1/(1._kp+A1*x))**2
    
  end function hf1i_epsilon_one


!epsilon_two(x)
  function hf1i_epsilon_two(x,A1)    
    implicit none
    real(kp) :: hf1i_epsilon_two
    real(kp), intent(in) :: x,A1
    
    hf1i_epsilon_two = 2._kp*hf1i_epsilon_one(x,A1) 
    
  end function hf1i_epsilon_two


!epsilon_three(x)
  function hf1i_epsilon_three(x,A1)    
    implicit none
    real(kp) :: hf1i_epsilon_three
    real(kp), intent(in) :: x,A1
    
    hf1i_epsilon_three = 2._kp*hf1i_epsilon_one(x,A1) 
    
  end function hf1i_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function hf1i_x_endinf(A1)
    implicit none
    real(kp), intent(in) :: A1
    real(kp) :: hf1i_x_endinf
    
    if (A1.le.0) then
    !Case when inflation proceeds from the left to the right
        hf1i_x_endinf = -sqrt(2._kp)-1._kp/A1
   
    else 
    !Case when inflation proceeds from the right to the left
    hf1i_x_endinf = sqrt(2._kp)-1._kp/A1

    end if

  end function hf1i_x_endinf



!this is integral[V(phi)/V'(phi) dphi]
  function hf1i_efold_primitive(x,A1)
    implicit none
    real(kp), intent(in) :: x,A1
    real(kp) :: hf1i_efold_primitive

    if (A1.eq.0._kp) stop 'hf1i_efold_primitive: A1=0!'

    hf1i_efold_primitive = 1._kp/(2._kp*A1) * (x+0.5_kp*A1*x**2)

  end function hf1i_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function hf1i_x_trajectory(bfold,xend,A1)
    implicit none
    real(kp), intent(in) :: bfold, A1, xend
    real(kp) :: hf1i_x_trajectory
    real(kp) :: OnePlusA1phiSquared

    OnePlusA1phiSquared=1._kp-2._kp*A1*(2._kp*A1*bfold-xend-A1/2._kp*xend**2)
	
    if (A1.le.0) then
    !Case when inflation proceeds from the left to the right
       hf1i_x_trajectory = (-1._kp-sqrt(OnePlusA1phiSquared))/A1
   
    else 
    !Case when inflation proceeds from the right to the left
       hf1i_x_trajectory = (-1._kp+sqrt(OnePlusA1phiSquared))/A1

    end if
    
           
  end function hf1i_x_trajectory


  
end module hf1isr
