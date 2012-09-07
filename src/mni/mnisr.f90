!slow-roll functions for the natural inflation with the minus sign potential
!
!V(phi) = M^4 [ 1 - cos(x/f) ]
!
!x = phi/Mp

module mnisr
  use infprec, only : kp
  implicit none

  private

  public  mni_norm_potential, mni_epsilon_one, mni_epsilon_two, mni_epsilon_three
  public  mni_x_endinf, mni_efold_primitive, mni_x_trajectory
  public  mni_norm_deriv_potential, mni_norm_deriv_second_potential
 
contains
!returns V/M^4
  function mni_norm_potential(x,f)
    implicit none
    real(kp) :: mni_norm_potential
    real(kp), intent(in) :: x,f

    mni_norm_potential = 1._kp-cos(x/f)

  end function mni_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function mni_norm_deriv_potential(x,f)
    implicit none
    real(kp) :: mni_norm_deriv_potential
    real(kp), intent(in) :: x,f

   mni_norm_deriv_potential = sin(x/f)/f

  end function mni_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function mni_norm_deriv_second_potential(x,f)
    implicit none
    real(kp) :: mni_norm_deriv_second_potential
    real(kp), intent(in) :: x,f

    mni_norm_deriv_second_potential =  cos(x/f)/(f**2)

  end function mni_norm_deriv_second_potential



!epsilon_one(x)
  function mni_epsilon_one(x,f)    
    implicit none
    real(kp) :: mni_epsilon_one
    real(kp), intent(in) :: x,f
    
    mni_epsilon_one = 1._kp/(2._kp*f**2) &
         *(sin(x/f))**2/(1._kp-cos(x/f))**2
    
  end function mni_epsilon_one


!epsilon_two(x)
  function mni_epsilon_two(x,f)
    implicit none
    real(kp) :: mni_epsilon_two
    real(kp), intent(in) :: x,f
    
    mni_epsilon_two = 2._kp/(f**2) &
         /(1._kp-cos(x/f))
    
  end function mni_epsilon_two


!epsilon_three(x)
  function mni_epsilon_three(x,f)
    implicit none
    real(kp) :: mni_epsilon_three
    real(kp), intent(in) :: x,f
    
    mni_epsilon_three = 2._kp*mni_epsilon_one(x,f) 
    
  end function mni_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function mni_x_endinf(f)
    implicit none
    real(kp), intent(in) :: f
    real(kp) :: mni_x_endinf
   
    mni_x_endinf = f*acos((-1._kp+2._kp*f**2)/(1._kp+2._kp*f**2))
   
  end function mni_x_endinf



!this is integral[V(phi)/V'(phi) dphi]
  function mni_efold_primitive(x,f)
    implicit none
    real(kp), intent(in) :: x,f
    real(kp) :: mni_efold_primitive
    
    mni_efold_primitive = -f**2 *2._kp &
         *log(abs(cos(x/(2._kp*f)) ))
         
   
  end function mni_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function mni_x_trajectory(bfold,xend,f)
    implicit none
    real(kp), intent(in) :: bfold,f,xend
    real(kp) :: mni_x_trajectory
    
    mni_x_trajectory = 2._kp*f*acos(cos(xend/(2._kp*f))*exp(bfold/(2._kp*f)))
       
  end function mni_x_trajectory

  
  
end module mnisr
