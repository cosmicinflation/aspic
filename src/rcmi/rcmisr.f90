!slow-roll functions for the radiatively corrected massive inflation potential
!
!V(phi) = M^4 x^2 [ 1 - 2 alpha x^2 ln(x) ]
!
!x = phi/Mp

module rcmisr
  use infprec, only : kp
  implicit none

  private

  public  rcmi_norm_potential, rcmi_epsilon_one, rcmi_epsilon_two
  public  rcmi_epsilon_three
  public  rcmi_x_endinf, rcmi_efold_primitive, rcmi_x_trajectory
 
contains
!returns V/M^4
  function rcmi_norm_potential(x,alpha)
    implicit none
    real(kp) :: rcmi_norm_potential
    real(kp), intent(in) :: x,alpha

    rcmi_norm_potential = x**2*(1._kp-2._kp*alpha*x**2*log(x))

  end function rcmi_norm_potential


!epsilon_one(x)
  function rcmi_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: rcmi_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    rcmi_epsilon_one = 2/(x**2) &
         *(1._kp-alpha*x**2-4._kp*alpha*x**2*log(x))**2 &
         /(1._kp-2._kp*alpha*x**2*log(x))**2
    
  end function rcmi_epsilon_one


!epsilon_two(x)
  function rcmi_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: rcmi_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    rcmi_epsilon_two = 4._kp/(x**2) &
         /(1._kp-2._kp*alpha*x**2*log(x))**2 &
         *(1._kp+3._kp*alpha*x**2 &
         -2._kp*alpha*x**2*log(x)+2._kp*alpha**2*x**4+2._kp*alpha**2*x**4*log(x) &
         +8._kp*alpha**2*x**4*(log(x))**2)
    
  end function rcmi_epsilon_two 


!epsilon_three(x)
  function rcmi_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: rcmi_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    rcmi_epsilon_three = 4._kp/(x**2) &
         /(1._kp-2._kp*alpha*x**2*log(x))**2 &
         *(1._kp-alpha*x**2-4._kp*alpha*x**2*log(x)) &
         *(1._kp-alpha*x**2-9._kp*alpha**2*x**4-4._kp*alpha**3*x**6 &
         -6._kp*alpha*x**2*log(x)-20._kp*alpha**2*x**4*log(x) &
         -6._kp*alpha**3*x**6*log(x) &
         -4._kp*alpha**3*x**6*(log(x))**2-16._kp*alpha**3*x**6*(log(x))**3) &
         *(1._kp+3._kp*alpha*x**2-2._kp*alpha*x**2*log(x)+2._kp*alpha**2*x**4 &
         +2._kp*alpha**2*x**4*log(x)+8._kp*alpha**2*x**4*(log(x))**2)**(-1)
    
  end function rcmi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function rcmi_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: rcmi_x_endinf
   
    rcmi_x_endinf = sqrt(2._kp)
   
  end function rcmi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function rcmi_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: rcmi_efold_primitive

    rcmi_efold_primitive = x**2/4._kp &
         +alpha**2/16._kp*x**4 &
         +alpha/4._kp*x**4*log(x)

  end function rcmi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function rcmi_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, xend, alpha
    real(kp) :: rcmi_x_trajectory
           
    rcmi_x_trajectory = 2._kp*sqrt(0.5_kp-bfold) &
         *(1._kp-alpha/2._kp*(0.5_kp-bfold) &
         -2._kp*alpha*(0.5_kp-bfold)*log(2._kp*sqrt(0.5_kp-bfold)) &
         +(1._kp+2._kp*log(2._kp))/(8._kp*(0.5_kp-bfold))*alpha)
    
   
  end function rcmi_x_trajectory

  
end module rcmisr
