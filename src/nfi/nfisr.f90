!slow-roll functions for the N-Formalism Inflation module
!
!V(phi) = M^4 exp(-a x^b)
!
!x = phi/Mp

module nficommon
  use infprec, only : kp
  implicit none

  private

  public  nfi_norm_potential, nfi_epsilon_one, nfi_epsilon_two
  public  nfi_epsilon_three
  public  nfi_x_endinf, nfi_efold_primitive, nfi_x_trajectory
  public  nfi_norm_deriv_potential, nfi_norm_deriv_second_potential
 
contains
!returns V/M^4
  function nfi_norm_potential(x,p)
    implicit none
    real(kp) :: nfi_norm_potential
    real(kp), intent(in) :: x,p

    nfi_norm_potential = x**p

  end function nfi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function nfi_norm_deriv_potential(x,p)
    implicit none
    real(kp) :: nfi_norm_deriv_potential
    real(kp), intent(in) :: x,p

   nfi_norm_deriv_potential = p*x**(p-1._kp)

  end function nfi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function nfi_norm_deriv_second_potential(x,p)
    implicit none
    real(kp) :: nfi_norm_deriv_second_potential
    real(kp), intent(in) :: x,p

    nfi_norm_deriv_second_potential =  p*(p-1._kp)*x**(p-2._kp)

  end function nfi_norm_deriv_second_potential


!epsilon_one(x)
  function nfi_epsilon_one(x,p)    
    implicit none
    real(kp) :: nfi_epsilon_one
    real(kp), intent(in) :: x,p
    
    nfi_epsilon_one = 0.5_kp*(p/x)**2
    
  end function nfi_epsilon_one


!epsilon_two(x)
  function nfi_epsilon_two(x,p)    
    implicit none
    real(kp) :: nfi_epsilon_two
    real(kp), intent(in) :: x,p
    
    nfi_epsilon_two = 2._kp*p/x**2
    
  end function nfi_epsilon_two


!epsilon_three(x)
  function nfi_epsilon_three(x,p)    
    implicit none
    real(kp) :: nfi_epsilon_three
    real(kp), intent(in) :: x,p
    
    nfi_epsilon_three = nfi_epsilon_two(x,p)
    
  end function nfi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function nfi_x_endinf(p)
    implicit none
    real(kp), intent(in) :: p
    real(kp) :: nfi_x_endinf
   
    nfi_x_endinf = p/sqrt(2._kp)
   
  end function nfi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function nfi_efold_primitive(x,p)
    implicit none
    real(kp), intent(in) :: x,p
    real(kp) :: nfi_efold_primitive

    if (p.eq.0._kp) stop 'nfi_efold_primitive: p=0!'

    nfi_efold_primitive = 0.5_kp*x**2/p

  end function nfi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function nfi_x_trajectory(bfold,xend,p)
    implicit none
    real(kp), intent(in) :: bfold, p, xend
    real(kp) :: nfi_x_trajectory
           
    nfi_x_trajectory = sqrt(-2._kp*p*bfold + xend**2)
    
   
  end function nfi_x_trajectory

  
end module nfisr
