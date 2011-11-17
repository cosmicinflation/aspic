!slow-roll functions for the large field potential
!
!V(phi) = M^4 (phi/Mp)^p
!
!x = phi/Mp

module lfisr
  use infprec, only : kp
  implicit none

  private

  public  lfi_norm_potential, lfi_epsilon_one, lfi_epsilon_two
  public  lfi_epsilon_three
  public  lfi_x_endinf, lfi_efold_primitive, lfi_x_trajectory
  public  lfi_norm_deriv_potential, lfi_norm_deriv_second_potential
 
contains
!returns V/M^4
  function lfi_norm_potential(x,p)
    implicit none
    real(kp) :: lfi_norm_potential
    real(kp), intent(in) :: x,p

    lfi_norm_potential = x**p

  end function lfi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function lfi_norm_deriv_potential(x,p)
    implicit none
    real(kp) :: lfi_norm_deriv_potential
    real(kp), intent(in) :: x,p

   lfi_norm_deriv_potential = p*x**(p-1._kp)

  end function lfi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function lfi_norm_deriv_second_potential(x,p)
    implicit none
    real(kp) :: lfi_norm_deriv_second_potential
    real(kp), intent(in) :: x,p

    lfi_norm_deriv_second_potential =  p*(p-1._kp)*x**(p-2._kp)

  end function lfi_norm_deriv_second_potential


!epsilon_one(x)
  function lfi_epsilon_one(x,p)    
    implicit none
    real(kp) :: lfi_epsilon_one
    real(kp), intent(in) :: x,p
    
    lfi_epsilon_one = 0.5_kp*(p/x)**2
    
  end function lfi_epsilon_one


!epsilon_two(x)
  function lfi_epsilon_two(x,p)    
    implicit none
    real(kp) :: lfi_epsilon_two
    real(kp), intent(in) :: x,p
    
    lfi_epsilon_two = 2._kp*p/x**2
    
  end function lfi_epsilon_two


!epsilon_three(x)
  function lfi_epsilon_three(x,p)    
    implicit none
    real(kp) :: lfi_epsilon_three
    real(kp), intent(in) :: x,p
    
    lfi_epsilon_three = lfi_epsilon_two(x,p)
    
  end function lfi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function lfi_x_endinf(p)
    implicit none
    real(kp), intent(in) :: p
    real(kp) :: lfi_x_endinf
   
    lfi_x_endinf = p/sqrt(2._kp)
   
  end function lfi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function lfi_efold_primitive(x,p)
    implicit none
    real(kp), intent(in) :: x,p
    real(kp) :: lfi_efold_primitive

    if (p.eq.0._kp) stop 'lfi_efold_primitive: p=0!'

    lfi_efold_primitive = 0.5_kp*x**2/p

  end function lfi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function lfi_x_trajectory(bfold,xend,p)
    implicit none
    real(kp), intent(in) :: bfold, p, xend
    real(kp) :: lfi_x_trajectory
           
    lfi_x_trajectory = sqrt(-2._kp*p*bfold + xend**2)
    
   
  end function lfi_x_trajectory

  
end module lfisr
