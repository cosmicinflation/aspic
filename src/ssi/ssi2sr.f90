!slow-roll functions for the sneutrino supersymmetric 2 potential
!
!
!V(phi) = M**4 [ 1 + alpha x**2 + beta x**4 ]
!
!2: alpha<0, beta<0
!
!x = phi/Mp

module ssi2sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ssicommon, only : ssi_norm_potential, ssi_norm_deriv_potential
  use ssicommon, only : ssi_norm_deriv_second_potential
  use ssicommon, only : ssi_epsilon_one, ssi_epsilon_two, ssi_epsilon_three
  use ssicommon, only : ssi_efold_primitive, find_ssitraj, ssi245_x_V_Equals_0


  implicit none

  private

  public ssi2_norm_potential
  public ssi2_epsilon_one, ssi2_epsilon_two, ssi2_epsilon_three
  public ssi2_x_endinf, ssi2_efold_primitive, ssi2_x_trajectory
  public ssi2_norm_deriv_potential, ssi2_norm_deriv_second_potential
  
contains

!returns V/M**4
  function ssi2_norm_potential(x,alpha,beta)    
    implicit none
    real(kp) :: ssi2_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssi2_norm_potential = ssi_norm_potential(x,alpha,beta)

  end function ssi2_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function ssi2_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi2_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta
    
    ssi2_norm_deriv_potential = ssi_norm_deriv_potential(x,alpha,beta)

  end function ssi2_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function ssi2_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi2_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta    

    ssi2_norm_deriv_second_potential &
         = ssi_norm_deriv_second_potential(x,alpha,beta)
    

  end function ssi2_norm_deriv_second_potential



!epsilon_one(x)
  function ssi2_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssi2_epsilon_one
    real(kp), intent(in) :: x,alpha,beta

    ssi2_epsilon_one = ssi_epsilon_one(x,alpha,beta)
    
  end function ssi2_epsilon_one


!epsilon_two(x)
  function ssi2_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssi2_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    ssi2_epsilon_two = ssi_epsilon_two(x,alpha,beta)
    
  end function ssi2_epsilon_two


!epsilon_three(x)
  function ssi2_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssi2_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
   
    ssi2_epsilon_three = ssi_epsilon_three(x,alpha,beta)
    
  end function ssi2_epsilon_three


!returns the position x where the potential vanishes
  function ssi2_x_V_Equals_0(alpha,beta)    
    implicit none
    real(kp) :: ssi2_x_V_Equals_0
    real(kp), intent(in) :: alpha,beta

    ssi2_x_V_Equals_0 = ssi245_x_V_Equals_0(alpha,beta)

    
  end function ssi2_x_V_Equals_0


!returns x at the end of inflation defined as epsilon1=1
  function ssi2_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: ssi2_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi2Data

    mini = epsilon(1._kp)
    maxi = ssi2_x_V_Equals_0(alpha,beta)*(1._kp-epsilon(1._kp))

!    print*,'ssi2_xend:  mini=',mini,'   maxi=',maxi,'   epsOne(mini)=',ssi2_epsilon_one(mini,alpha,beta), &
!                 '   epsOne(maxi)=',ssi2_epsilon_one(maxi,alpha,beta)
!    pause

    ssi2Data%real1 = alpha
    ssi2Data%real2 = beta
    
    ssi2_x_endinf = zbrent(find_ssi2endinf,mini,maxi,tolFind,ssi2Data)

  end function ssi2_x_endinf



  function find_ssi2endinf(x,ssi2Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssi2Data
    real(kp) :: find_ssi2endinf
    real(kp) :: alpha,beta

    alpha = ssi2Data%real1
    beta = ssi2Data%real2
    
    find_ssi2endinf = ssi2_epsilon_one(x,alpha,beta) - 1._kp
   
  end function find_ssi2endinf

!this is integral[V(phi)/V'(phi) dphi]
  function ssi2_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssi2_efold_primitive

    ssi2_efold_primitive = ssi_efold_primitive(x,alpha,beta)
  
  end function ssi2_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ssi2_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend,beta
    real(kp) :: ssi2_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi2Data


    mini = epsilon(1._kp)
    maxi = ssi2_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))
  
    ssi2Data%real1 = alpha
    ssi2Data%real2 = beta
    ssi2Data%real3 = -bfold + ssi2_efold_primitive(xend,alpha,beta)
    
    ssi2_x_trajectory = zbrent(find_ssitraj,mini,maxi,tolFind,ssi2Data)
       
  end function ssi2_x_trajectory
 

end module ssi2sr
