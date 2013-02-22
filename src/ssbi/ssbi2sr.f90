!slow-roll functions for the spontaneous symmetry breaking 2 potential
!
!
!V(phi) = M**4 [ 1 + alpha x**2 + beta x**4 ]
!
!2: alpha<0, beta<0
!
!x = phi/Mp

module ssbi2sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ssbicommon, only : ssbi_norm_potential, ssbi_norm_deriv_potential
  use ssbicommon, only : ssbi_norm_deriv_second_potential
  use ssbicommon, only : ssbi_epsilon_one, ssbi_epsilon_two, ssbi_epsilon_three
  use ssbicommon, only : ssbi_efold_primitive, find_ssbi_x_trajectory, ssbi245_x_potzero


  implicit none

  private

  public ssbi2_norm_potential
  public ssbi2_epsilon_one, ssbi2_epsilon_two, ssbi2_epsilon_three
  public ssbi2_x_endinf, ssbi2_efold_primitive, ssbi2_x_trajectory
  public ssbi2_norm_deriv_potential, ssbi2_norm_deriv_second_potential
  
contains

!returns V/M**4
  function ssbi2_norm_potential(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi2_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssbi2_norm_potential = ssbi_norm_potential(x,alpha,beta)

  end function ssbi2_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function ssbi2_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi2_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta
    
    ssbi2_norm_deriv_potential = ssbi_norm_deriv_potential(x,alpha,beta)

  end function ssbi2_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function ssbi2_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi2_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta    

    ssbi2_norm_deriv_second_potential &
         = ssbi_norm_deriv_second_potential(x,alpha,beta)
    

  end function ssbi2_norm_deriv_second_potential



!epsilon_one(x)
  function ssbi2_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi2_epsilon_one
    real(kp), intent(in) :: x,alpha,beta

    ssbi2_epsilon_one = ssbi_epsilon_one(x,alpha,beta)
    
  end function ssbi2_epsilon_one


!epsilon_two(x)
  function ssbi2_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi2_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    ssbi2_epsilon_two = ssbi_epsilon_two(x,alpha,beta)
    
  end function ssbi2_epsilon_two


!epsilon_three(x)
  function ssbi2_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi2_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
   
    ssbi2_epsilon_three = ssbi_epsilon_three(x,alpha,beta)
    
  end function ssbi2_epsilon_three


!returns the position x where the potential vanishes
  function ssbi2_x_potzero(alpha,beta)    
    implicit none
    real(kp) :: ssbi2_x_potzero
    real(kp), intent(in) :: alpha,beta

    ssbi2_x_potzero = ssbi245_x_potzero(alpha,beta)

  end function ssbi2_x_potzero


!returns x at the end of inflation defined as epsilon1=1
  function ssbi2_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: ssbi2_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi2Data

    mini = epsilon(1._kp)
    maxi = ssbi2_x_potzero(alpha,beta)*(1._kp-epsilon(1._kp))

!    print*,'ssbi2_xend:  mini=',mini,'   maxi=',maxi,'   epsOne(mini)=',ssbi2_epsilon_one(mini,alpha,beta), &
!                 '   epsOne(maxi)=',ssbi2_epsilon_one(maxi,alpha,beta)
!    pause

    ssbi2Data%real1 = alpha
    ssbi2Data%real2 = beta
    
    ssbi2_x_endinf = zbrent(find_ssbi2_x_endinf,mini,maxi,tolFind,ssbi2Data)

  end function ssbi2_x_endinf



  function find_ssbi2_x_endinf(x,ssbi2Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssbi2Data
    real(kp) :: find_ssbi2_x_endinf
    real(kp) :: alpha,beta

    alpha = ssbi2Data%real1
    beta = ssbi2Data%real2
    
    find_ssbi2_x_endinf = ssbi2_epsilon_one(x,alpha,beta) - 1._kp
   
  end function find_ssbi2_x_endinf

!this is integral[V(phi)/V'(phi) dphi]
  function ssbi2_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssbi2_efold_primitive

    ssbi2_efold_primitive = ssbi_efold_primitive(x,alpha,beta)
  
  end function ssbi2_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ssbi2_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend,beta
    real(kp) :: ssbi2_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi2Data


    mini = epsilon(1._kp)
    maxi = ssbi2_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))
  
    ssbi2Data%real1 = alpha
    ssbi2Data%real2 = beta
    ssbi2Data%real3 = -bfold + ssbi2_efold_primitive(xend,alpha,beta)
    
    ssbi2_x_trajectory = zbrent(find_ssbi_x_trajectory,mini,maxi,tolFind,ssbi2Data)
       
  end function ssbi2_x_trajectory
 

end module ssbi2sr
