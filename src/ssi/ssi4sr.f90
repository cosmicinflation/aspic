!slow-roll functions for the sneutrino supersymmetric 4 potential
!
!
!V(phi) = M**4 [ 1 + alpha x**2 + beta x**4 ]
!
!2: alpha>0, beta<0m x^2 > -alpha / ( 2 beta )
!
!x = phi/Mp

module ssi4sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ssicommon, only : ssi_norm_potential, ssi_norm_deriv_potential, ssi_norm_deriv_second_potential
  use ssicommon, only : ssi_epsilon_one, ssi_epsilon_two, ssi_epsilon_three
  use ssicommon, only : ssi_efold_primitive, find_ssitraj, ssi245_x_V_Equals_0, ssi3456_x_Vprime_Equals_0


  implicit none

  private

  public ssi4_norm_potential
  public ssi4_epsilon_one, ssi4_epsilon_two, ssi4_epsilon_three
  public ssi4_x_endinf, ssi4_efold_primitive, ssi4_x_trajectory
  public ssi4_norm_deriv_potential, ssi4_norm_deriv_second_potential
  
contains

!returns V/M**4
  function ssi4_norm_potential(x,alpha,beta)    
    implicit none
    real(kp) :: ssi4_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssi4_norm_potential = ssi_norm_potential(x,alpha,beta)

  end function ssi4_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function ssi4_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi4_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta
    
    ssi4_norm_deriv_potential = ssi_norm_deriv_potential(x,alpha,beta)

  end function ssi4_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function ssi4_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi4_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta    

    ssi4_norm_deriv_second_potential &
         = ssi_norm_deriv_second_potential(x,alpha,beta)
    

  end function ssi4_norm_deriv_second_potential



!epsilon_one(x)
  function ssi4_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssi4_epsilon_one
    real(kp), intent(in) :: x,alpha,beta

    ssi4_epsilon_one = ssi_epsilon_one(x,alpha,beta)
    
  end function ssi4_epsilon_one


!epsilon_two(x)
  function ssi4_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssi4_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    ssi4_epsilon_two = ssi_epsilon_two(x,alpha,beta)
    
  end function ssi4_epsilon_two


!epsilon_three(x)
  function ssi4_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssi4_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
   
    ssi4_epsilon_three = ssi_epsilon_three(x,alpha,beta)
    
  end function ssi4_epsilon_three


!returns the position x where the potential vanishes
  function ssi4_x_V_Equals_0(alpha,beta)    
    implicit none
    real(kp) :: ssi4_x_V_Equals_0
    real(kp), intent(in) :: alpha,beta

    ssi4_x_V_Equals_0 = ssi245_x_V_Equals_0(alpha,beta)

    
  end function ssi4_x_V_Equals_0


!returns x at the end of inflation defined as epsilon1=1
  function ssi4_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: ssi4_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi4Data

    mini = ssi3456_x_Vprime_Equals_0(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssi4_x_V_Equals_0(alpha,beta)*(1._kp-epsilon(1._kp))

    ssi4Data%real1 = alpha
    ssi4Data%real2 = beta
    
    ssi4_x_endinf = zbrent(find_ssi4endinf,mini,maxi,tolFind,ssi4Data)

  end function ssi4_x_endinf



  function find_ssi4endinf(x,ssi4Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssi4Data
    real(kp) :: find_ssi4endinf
    real(kp) :: alpha,beta

    alpha = ssi4Data%real1
    beta = ssi4Data%real2
    
    find_ssi4endinf = ssi4_epsilon_one(x,alpha,beta) - 1._kp
   
  end function find_ssi4endinf

!this is integral[V(phi)/V'(phi) dphi]
  function ssi4_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssi4_efold_primitive

    ssi4_efold_primitive = ssi_efold_primitive(x,alpha,beta)
  
  end function ssi4_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ssi4_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend,beta
    real(kp) :: ssi4_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi4Data


    mini = ssi3456_x_Vprime_Equals_0(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssi4_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))
  
    ssi4Data%real1 = alpha
    ssi4Data%real2 = beta
    ssi4Data%real3 = -bfold + ssi4_efold_primitive(xend,alpha,beta)
    
    ssi4_x_trajectory = zbrent(find_ssitraj,mini,maxi,tolFind,ssi4Data)
       
  end function ssi4_x_trajectory
 

end module ssi4sr
