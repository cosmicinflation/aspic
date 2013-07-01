!slow-roll functions for the spontaneous symmetry breaking 4 potential
!
!
!V(phi) = M**4 [ 1 + alpha x**2 + beta x**4 ]
!
!2: alpha>0, beta<0  x^2 > -alpha / ( 2 beta )
!
!x = phi/Mp

module ssbi4sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ssbicommon, only : ssbi_norm_potential, ssbi_norm_deriv_potential
  use ssbicommon, only : ssbi_norm_deriv_second_potential, ssbi3456_x_derivpotzero
  use ssbicommon, only : ssbi_epsilon_one, ssbi_epsilon_two, ssbi_epsilon_three
  use ssbicommon, only : ssbi_efold_primitive, find_ssbi_x_trajectory, ssbi245_x_potzero


  implicit none

  private

  public ssbi4_norm_potential
  public ssbi4_epsilon_one, ssbi4_epsilon_two, ssbi4_epsilon_three
  public ssbi4_x_endinf, ssbi4_efold_primitive, ssbi4_x_trajectory
  public ssbi4_norm_deriv_potential, ssbi4_norm_deriv_second_potential
  public ssbi4_x_potmax, ssbi4_x_potzero, ssbi4_eps2_x_potmax
  
contains

!returns V/M**4
  function ssbi4_norm_potential(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi4_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssbi4_norm_potential = ssbi_norm_potential(x,alpha,beta)

  end function ssbi4_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function ssbi4_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi4_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta
    
    ssbi4_norm_deriv_potential = ssbi_norm_deriv_potential(x,alpha,beta)

  end function ssbi4_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function ssbi4_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi4_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta    

    ssbi4_norm_deriv_second_potential &
         = ssbi_norm_deriv_second_potential(x,alpha,beta)
    

  end function ssbi4_norm_deriv_second_potential



!epsilon_one(x)
  function ssbi4_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi4_epsilon_one
    real(kp), intent(in) :: x,alpha,beta

    ssbi4_epsilon_one = ssbi_epsilon_one(x,alpha,beta)
    
  end function ssbi4_epsilon_one


!epsilon_two(x)
  function ssbi4_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi4_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    ssbi4_epsilon_two = ssbi_epsilon_two(x,alpha,beta)
    
  end function ssbi4_epsilon_two


!epsilon_three(x)
  function ssbi4_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi4_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
   
    ssbi4_epsilon_three = ssbi_epsilon_three(x,alpha,beta)
    
  end function ssbi4_epsilon_three


!returns the position x where the potential vanishes
  function ssbi4_x_potzero(alpha,beta)    
    implicit none
    real(kp) :: ssbi4_x_potzero
    real(kp), intent(in) :: alpha,beta

    ssbi4_x_potzero = ssbi245_x_potzero(alpha,beta)

    
  end function ssbi4_x_potzero


!returns the position x where the potential is maximal
  function ssbi4_x_potmax(alpha,beta)    
    implicit none
    real(kp) :: ssbi4_x_potmax
    real(kp), intent(in) :: alpha,beta

    ssbi4_x_potmax = ssbi3456_x_derivpotzero(alpha,beta)

    
  end function ssbi4_x_potmax

! returns the value of eps2 at the maximum of the potential
  function ssbi4_eps2_x_potmax(alpha,beta)    
    implicit none
    real(kp) :: ssbi4_eps2_x_potmax
    real(kp), intent(in) :: alpha,beta

    ssbi4_eps2_x_potmax = ssbi4_epsilon_two(ssbi4_x_potmax(alpha,beta),alpha,beta)
    
  end function ssbi4_eps2_x_potmax


!returns x at the end of inflation defined as epsilon1=1
  function ssbi4_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: ssbi4_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi4Data

    mini = ssbi4_x_potmax(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssbi4_x_potzero(alpha,beta)*(1._kp-epsilon(1._kp))

    ssbi4Data%real1 = alpha
    ssbi4Data%real2 = beta
    
    ssbi4_x_endinf = zbrent(find_ssbi4_x_endinf,mini,maxi,tolFind,ssbi4Data)

  end function ssbi4_x_endinf



  function find_ssbi4_x_endinf(x,ssbi4Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssbi4Data
    real(kp) :: find_ssbi4_x_endinf
    real(kp) :: alpha,beta

    alpha = ssbi4Data%real1
    beta = ssbi4Data%real2
    
    find_ssbi4_x_endinf = ssbi4_epsilon_one(x,alpha,beta) - 1._kp
   
  end function find_ssbi4_x_endinf

!this is integral[V(phi)/V'(phi) dphi]
  function ssbi4_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssbi4_efold_primitive

    ssbi4_efold_primitive = ssbi_efold_primitive(x,alpha,beta)
  
  end function ssbi4_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ssbi4_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend,beta
    real(kp) :: ssbi4_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi4Data


    mini = ssbi3456_x_derivpotzero(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssbi4_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))
  
    ssbi4Data%real1 = alpha
    ssbi4Data%real2 = beta
    ssbi4Data%real3 = -bfold + ssbi4_efold_primitive(xend,alpha,beta)
    
    ssbi4_x_trajectory = zbrent(find_ssbi_x_trajectory,mini,maxi,tolFind,ssbi4Data)
       
  end function ssbi4_x_trajectory
 

end module ssbi4sr
