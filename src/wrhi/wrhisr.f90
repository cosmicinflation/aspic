!slow-roll functions for the Witten-O'Raifeartaigh Higgs inflation potential
!
!V(phi) = M^4 Ln[x]^2
!
!x=phi/phi0
!phi0=phi0/Mp


module wrhisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  use inftools, only : zbrent
  implicit none

  private

  public  wrhi_norm_potential, wrhi_epsilon_one, wrhi_epsilon_two, wrhi_epsilon_three
  public  wrhi_x_endinf, wrhi_efold_primitive, wrhi_x_trajectory
  public  wrhi_norm_deriv_potential, wrhi_norm_deriv_second_potential
 
contains
!returns V/M**4
  function wrhi_norm_potential(x,phi0)
    implicit none
    real(kp) :: wrhi_norm_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: phi0

    wrhi_norm_potential = log(x)**2

  end function wrhi_norm_potential



!returns the first derivative of the potential with respect to x=phi/phi0, divided by M**4
  function wrhi_norm_deriv_potential(x,phi0)
    implicit none
    real(kp) :: wrhi_norm_deriv_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: phi0

   wrhi_norm_deriv_potential = 2._kp*log(x)/x

  end function wrhi_norm_deriv_potential



!returns the second derivative of the potential with respect to x=phi/phi0, divided by M**4
  function wrhi_norm_deriv_second_potential(x,phi0)
    implicit none
    real(kp) :: wrhi_norm_deriv_second_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: phi0

    wrhi_norm_deriv_second_potential = 2._kp*(1._kp-log(x))/(x**2)

  end function wrhi_norm_deriv_second_potential



!epsilon_one(x)
  function wrhi_epsilon_one(x,phi0)    
    implicit none
    real(kp) :: wrhi_epsilon_one
    real(kp), intent(in) :: x,phi0

  
    wrhi_epsilon_one = 2._kp/((phi0*x*log(x))**2)

    
  end function wrhi_epsilon_one


!epsilon_two(x)
  function wrhi_epsilon_two(x,phi0)    
    implicit none
    real(kp) :: wrhi_epsilon_two
    real(kp), intent(in) :: x,phi0
    
    wrhi_epsilon_two = 4._kp*(1._kp+log(x))/((phi0*x*log(x))**2)
 
  end function wrhi_epsilon_two


!epsilon_three(x)
  function wrhi_epsilon_three(x,phi0)    
    implicit none
    real(kp) :: wrhi_epsilon_three
    real(kp), intent(in) :: x,phi0
    
    wrhi_epsilon_three = 2._kp*(2._kp+3._kp*log(x)+2._kp*log(x)**2)/ &
                         ((phi0*x*log(x))**2*(1._kp+log(x)))
   
  end function wrhi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function wrhi_x_endinf(phi0)
    implicit none
    real(kp), intent(in) :: phi0
    real(kp) :: wrhi_x_endinf

    wrhi_x_endinf = exp(lambert(sqrt(2._kp)/phi0,0))


  end function wrhi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function wrhi_efold_primitive(x,phi0)
    implicit none
    real(kp), intent(in) :: x,phi0
    real(kp) :: wrhi_efold_primitive


    wrhi_efold_primitive = 0.25_kp*phi0**2*x**2*(log(x)-0.5_kp)


  end function wrhi_efold_primitive

!returns x=phi/mi at bfold=-efolds before the end of inflation, ie N-Nend
  function wrhi_x_trajectory(bfold,xend,phi0)
    implicit none
    real(kp), intent(in) :: bfold, xend, phi0
    real(kp)  :: wrhi_x_trajectory

    wrhi_x_trajectory = exp(0.5*lambert(8._kp/(exp(1.)*phi0**2)*(-bfold)+ &
                        2._kp/exp(1._kp)*xend**2*log(xend)- &
                        xend**2/exp(1._kp),0)+0.5_kp)
       
  end function wrhi_x_trajectory


end module wrhisr
