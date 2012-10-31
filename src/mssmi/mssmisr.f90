!slow-roll functions for the MSSM inflation potential
!
!V(phi) = M^4 [ x^2 - alpha x^6 + 9/20 alpha^2 x^10 ]
!
!x = phi/Mp

module mssmisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : atanh
  use mssmicommon, only : gmssmi_gen_norm_potential, gmssmi_gen_norm_deriv_potential, gmssmi_gen_norm_deriv_second_potential
  use mssmicommon, only : gmssmi_gen_epsilon_one, gmssmi_gen_epsilon_two, gmssmi_gen_epsilon_three
  use mssmicommon, only : gmssmi_gen_x_endinf, gmssmi_gen_x_epsilon1_min
  implicit none

  private

  public  mssmi_norm_potential, mssmi_epsilon_one, mssmi_epsilon_two, mssmi_epsilon_three
  public  mssmi_x_endinf, mssmi_efold_primitive, mssmi_x_trajectory
  public  mssmi_norm_deriv_potential, mssmi_norm_deriv_second_potential
  public  mssmi_x_epsilon1_min

 
contains


  function mssmi_beta(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: mssmi_beta

    mssmi_beta = 9._kp/20._kp*alpha**2.

  end function mssmi_beta


!returns V/M**4
  function mssmi_norm_potential(x,alpha)
    implicit none
    real(kp) :: mssmi_norm_potential
    real(kp), intent(in) :: x,alpha

    mssmi_norm_potential = gmssmi_gen_norm_potential(x,alpha,mssmi_beta(alpha))

  end function mssmi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function mssmi_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: mssmi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

   mssmi_norm_deriv_potential = gmssmi_gen_norm_deriv_potential(x,alpha,mssmi_beta(alpha))

  end function mssmi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function mssmi_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: mssmi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    mssmi_norm_deriv_second_potential = gmssmi_gen_norm_deriv_second_potential(x,alpha,mssmi_beta(alpha))

  end function mssmi_norm_deriv_second_potential



!epsilon_one(x)
  function mssmi_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: mssmi_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    mssmi_epsilon_one = gmssmi_gen_epsilon_one(x,alpha,mssmi_beta(alpha))
    
  end function mssmi_epsilon_one


!epsilon_two(x)
  function mssmi_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: mssmi_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    mssmi_epsilon_two = gmssmi_gen_epsilon_two(x,alpha,mssmi_beta(alpha)) 
    
  end function mssmi_epsilon_two


!epsilon_three(x)
  function mssmi_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: mssmi_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    mssmi_epsilon_three = gmssmi_gen_epsilon_three(x,alpha,mssmi_beta(alpha)) 
    
  end function mssmi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function mssmi_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: mssmi_x_endinf
    
    mssmi_x_endinf = gmssmi_gen_x_endinf(alpha,mssmi_beta(alpha))
   

  end function mssmi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function mssmi_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: mssmi_efold_primitive

    mssmi_efold_primitive = x**2/20._kp+2._kp/15._kp*x**2/(2._kp-3._kp*alpha*x**4) &
               +2._kp/15._kp*sqrt(2._kp/(3._kp*alpha))*real(atanh(sqrt(3._kp*alpha/2._kp)*x**2))

  end function mssmi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function mssmi_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: mssmi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: mssmiData

  
    mini = xend

    maxi = mssmi_x_epsilon1_min(alpha)*10._kp !Scales as the position of the inflexion point

    mssmiData%real1 = alpha
    mssmiData%real2 = -bfold + mssmi_efold_primitive(xend,alpha)
    
    mssmi_x_trajectory = zbrent(find_mssmitraj,mini,maxi,tolFind,mssmiData)
       
  end function mssmi_x_trajectory

  function find_mssmitraj(x,mssmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: mssmiData
    real(kp) :: find_mssmitraj
    real(kp) :: alpha,NplusNuend

    alpha= mssmiData%real1
    NplusNuend = mssmiData%real2

    find_mssmitraj = mssmi_efold_primitive(x,alpha) - NplusNuend
   
  end function find_mssmitraj



!Returns the position of the first local minimum of epsilon1
  function mssmi_x_epsilon1_min(alpha)   
    implicit none
    real(kp) :: mssmi_x_epsilon1_min
    real(kp), intent(in) :: alpha

	mssmi_x_epsilon1_min = gmssmi_gen_x_epsilon1_min(alpha,mssmi_beta(alpha)) 
    
  end function mssmi_x_epsilon1_min

  
end module mssmisr
