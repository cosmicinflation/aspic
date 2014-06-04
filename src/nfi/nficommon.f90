!common slow-roll functions for the N-Formalism Inflation module
!
!V(phi) = M^4 exp(-a x^b)
!
!x = phi/Mp

module nficommon
  use infprec, only : kp
  implicit none

  private

  public  nfi_norm_potential
  public  nfi_norm_deriv_potential, nfi_norm_deriv_second_potential
  public  nfi_epsilon_one, nfi_epsilon_two, nfi_epsilon_three
  public  nfi_x_epsoneunity, nfi_efold_primitive, nfi_x_trajectory
  
 
contains
!returns V/M^4
  function nfi_norm_potential(x,a,b)
    implicit none
    real(kp) :: nfi_norm_potential
    real(kp), intent(in) :: x,a,b

    nfi_norm_potential = exp(-a*x**b)

  end function nfi_norm_potential


!first derivative of the potential/M^4 with respect to x
  function nfi_norm_deriv_potential(x,a,b)
    implicit none
    real(kp) :: nfi_norm_deriv_potential
    real(kp), intent(in) :: x,a,b

   nfi_norm_deriv_potential = -a*b*x**(b-1._kp)*exp(-a*x**b)

  end function nfi_norm_deriv_potential


!second derivative of the potential/M^4 with respect to x
  function nfi_norm_deriv_second_potential(x,a,b)
    implicit none
    real(kp) :: nfi_norm_deriv_second_potential
    real(kp), intent(in) :: x,a,b

    nfi_norm_deriv_second_potential =  (1._kp-b + a*b*x**b) &
         *a*b*x**(b-2._kp)*exp(-a*x**b)

  end function nfi_norm_deriv_second_potential


!epsilon_one(x)
  function nfi_epsilon_one(x,a,b)    
    implicit none
    real(kp) :: nfi_epsilon_one
    real(kp), intent(in) :: x,a,b
    
    nfi_epsilon_one = 0.5_kp*a*a*b*b*(x*x)**(b-1._kp)
    
  end function nfi_epsilon_one


!epsilon_two(x)
  function nfi_epsilon_two(x,a,b)    
    implicit none
    real(kp) :: nfi_epsilon_two
    real(kp), intent(in) :: x,a,b
    
    nfi_epsilon_two = 2._kp*a*b*(b - 1._kp)*x**(b-2._kp)
    
  end function nfi_epsilon_two


!epsilon_three(x)
  function nfi_epsilon_three(x,a,b)    
    implicit none
    real(kp) :: nfi_epsilon_three
    real(kp), intent(in) :: x,a,b
    
    nfi_epsilon_three = a*b*(b-2._kp)*x**(b-2._kp)
    
  end function nfi_epsilon_three


!returns a vector containing the two roots of eps1=1
  function nfi_x_epsoneunity(a,b)
    implicit none
    real(kp), intent(in) :: a,b
    real(kp), dimension(2) :: nfi_x_epsoneunity
   
    if (b.eq.1._kp) then
       stop 'nfi_x_epsoneunity: b=1 is PLI'
    elseif (b.eq.0._kp) then
       stop 'nfi_x_epsoneunity: b=0 is singular'
    endif

!positive root
    nfi_x_epsoneunity(1) = (2._kp/a/a/b/b)**(0.5_kp/(b-1._kp))
!negative root
    nfi_x_epsoneunity(2) = - nfi_x_epsoneunity(1)
   
  end function nfi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function nfi_efold_primitive(x,a,b)
    implicit none
    real(kp), intent(in) :: x,a,b
    real(kp) :: nfi_efold_primitive

    if (b.eq.2._kp) then
       nfi_efold_primitive = -0.5_kp*log(x)/a
    else
       nfi_efold_primitive = x**(2._kp-b)/(a*b*(b-2._kp))
    endif

  end function nfi_efold_primitive


!Return x as a function of bfold = N - N_end
  function nfi_x_trajectory(bfold,a,b,xend)
    implicit none
    real(kp) :: nfi_x_trajectory
    real(kp), intent(in) :: bfold,a,b,xend

    if (b.eq.2._kp) then
       nfi_x_trajectory = xend*exp(2._kp*a*bfold)
    else 
       nfi_x_trajectory = (xend**(2._kp-b) + a*b*(2._kp-b)*bfold) &
            **(1._kp/(2._kp-b))
    endif

  end function nfi_x_trajectory


end module nficommon
