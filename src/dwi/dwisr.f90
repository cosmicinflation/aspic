!slow-roll functions for the double well inflation potential
!
!V(phi) = M^4 * [(phi/pi0)^2 - 1]^2
!
!x=phi/phi0
!phi0=phi0/Mp


module dwisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  implicit none

  private

  public  dwi_norm_potential, dwi_epsilon_one, dwi_epsilon_two, dwi_epsilon_three
  public  dwi_x_endinf, dwi_efold_primitive, dwi_x_trajectory
  public  dwi_norm_deriv_potential, dwi_norm_deriv_second_potential
 
contains
!returns V/M^4
  function dwi_norm_potential(x,phi0)
    implicit none
    real(kp) :: dwi_norm_potential
    real(kp), intent(in) :: x,phi0

    dwi_norm_potential = (x**2-1._kp)**2

  end function dwi_norm_potential



!returns the first derivative of the potential with respect to phi, divided by M^4
  function dwi_norm_deriv_potential(x,phi0)
    implicit none
    real(kp) :: dwi_norm_deriv_potential
    real(kp), intent(in) :: x,phi0

   dwi_norm_deriv_potential = (4._kp*x*(-1._kp+x**2))/phi0

  end function dwi_norm_deriv_potential



!returns the second derivative of the potential with respect to phi, divided by M^4
  function dwi_norm_deriv_second_potential(x,phi0)
    implicit none
    real(kp) :: dwi_norm_deriv_second_potential
    real(kp), intent(in) :: x,phi0

    dwi_norm_deriv_second_potential = (4._kp*(-1._kp+3._kp*x**2))/phi0**2

  end function dwi_norm_deriv_second_potential



!epsilon_one(x)
  function dwi_epsilon_one(x,phi0)    
    implicit none
    real(kp) :: dwi_epsilon_one
    real(kp), intent(in) :: x,phi0

  
    dwi_epsilon_one = (8._kp*x**2)/(phi0**2*(-1._kp+x**2)**2)

    
  end function dwi_epsilon_one


!epsilon_two(x)
  function dwi_epsilon_two(x,phi0)    
    implicit none
    real(kp) :: dwi_epsilon_two
    real(kp), intent(in) :: x,phi0
    
    dwi_epsilon_two = (8._kp*(1._kp+x**2))/(phi0**2*(-1._kp+x**2)**2)
    
  end function dwi_epsilon_two


!epsilon_three(x)
  function dwi_epsilon_three(x,phi0)    
    implicit none
    real(kp) :: dwi_epsilon_three
    real(kp), intent(in) :: x,phi0
    
    dwi_epsilon_three = (8._kp*x**2*(3._kp+x**2))/(phi0**2*(-1._kp+x**2)**2*(1._kp+x**2))
   
  end function dwi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function dwi_x_endinf(phi0)
    implicit none
    real(kp), intent(in) :: phi0
    real(kp) :: dwi_x_endinf
    

    dwi_x_endinf = -sqrt(2._kp)/phi0+sqrt(1._kp+2._kp/phi0**2)


  end function dwi_x_endinf




!this is integral[V(phi)/V'(phi) dphi]
  function dwi_efold_primitive(x,phi0)
    implicit none
    real(kp), intent(in) :: x,phi0
    real(kp) :: dwi_efold_primitive


        dwi_efold_primitive = 0.25_kp*phi0**2*(0.5_kp*x**2-log(x))


  end function dwi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function dwi_x_trajectory(bfold,xend,phi0)
    implicit none
    real(kp), intent(in) :: bfold, phi0, xend
    real(kp) :: dwi_x_trajectory


    dwi_x_trajectory = sqrt(-lambert(-xend**2*exp(-xend**2+8._kp*bfold/phi0**2),0))

       
  end function dwi_x_trajectory


end module dwisr
