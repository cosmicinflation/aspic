!slow-roll function for string axion II inflation at decreasing field
!values x, when the potential is an always monotonous function of x,
!i.e., does not have any minimum or maximum.
!
!V(phi) = M^4 [1 - cos(x) + alpha x sin(x) + (1/2) alpha beta x^2]
!
!x=phi/mu
!
!
!
module saiii3sr
  use infprec, only : kp, pi, tolkp, transfert
  use inftools, only : zbrent
  use saiiicommon, only : saiii_norm_potential, saiii_norm_deriv_potential
  use saiiicommon, only : saiii_norm_deriv_second_potential
  use saiiicommon, only : saiii_epsilon_one, saiii_epsilon_two, saiii_epsilon_three
  use saiiicommon, only : saiii_x_potzero, saiii_x_potmax, saiii_efold_primitive
  use saiiicommon, only : saiii_x_epsoneunity, find_saiii_x_trajectory
  use saiiicommon, only : saiii_check_minima, saiiiXBig

  
  implicit none


  
  private

  public saiii3_norm_potential, saiii3_norm_deriv_potential, saiii3_norm_deriv_second_potential
  public saiii3_epsilon_one, saiii3_epsilon_two, saiii3_epsilon_three
  public saiii3_x_endinf, saiii3_x_trajectory, saiii3_efold_primitive
  
contains

!check params for no extremum
  function saiii3_check_params(alpha,beta,mu)
    implicit none
    logical :: saiii3_check_params
    real(kp), intent(in) :: alpha,beta,mu

    saiii3_check_params = .not.saiii_check_minima(alpha,beta,mu)
    
  end function saiii3_check_params

  
  function saiii3_norm_potential(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii3_norm_potential
    real(kp), intent(in) :: x,alpha,beta,mu

    saiii3_norm_potential = saiii_norm_potential(x,alpha,beta,mu)
    
  end function saiii3_norm_potential


  
!derivative with respect to x (not phi!)  
  function saiii3_norm_deriv_potential(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii3_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta,mu

    saiii3_norm_deriv_potential =  saiii_norm_deriv_potential(x,alpha,beta,mu)

  end function saiii3_norm_deriv_potential

  

!second derivative with respect to x (not phi!)    
  function saiii3_norm_deriv_second_potential(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii3_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta,mu
    
    saiii3_norm_deriv_second_potential = saiii_norm_deriv_second_potential(x,alpha,beta,mu)

  end function saiii3_norm_deriv_second_potential



  
  function saiii3_epsilon_one(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii3_epsilon_one
    real(kp), intent(in) :: x,alpha,beta,mu
    
    saiii3_epsilon_one = saiii_epsilon_one(x,alpha,beta,mu)
    
  end function saiii3_epsilon_one
 
  

  
  function saiii3_epsilon_two(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii3_epsilon_two
    real(kp), intent(in) :: x,alpha,beta,mu
    
    saiii3_epsilon_two = saiii_epsilon_two(x,alpha,beta,mu)
    
  end function saiii3_epsilon_two



  
  function saiii3_epsilon_three(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii3_epsilon_three
    real(kp), intent(in) :: x,alpha,beta,mu
    
    saiii3_epsilon_three = saiii_epsilon_three(x,alpha,beta,mu)

  end function saiii3_epsilon_three



  
  function saiii3_x_endinf(alpha,beta,mu)
    implicit none
    real(kp) :: saiii3_x_endinf
    real(kp), intent(in) :: alpha, beta, mu

    real(kp), dimension(2) :: xepsone
    
    xepsone = saiii_x_epsoneunity(alpha,beta,mu)
    
    saiii3_x_endinf = xepsone(1)
    
  end function saiii3_x_endinf


  
  
  function saiii3_efold_primitive(x,alpha,beta,mu)
    implicit none
    real(kp), intent(in) :: x, alpha,beta,mu
    real(kp) :: saiii3_efold_primitive   
    
    saiii3_efold_primitive = saiii_efold_primitive(x,alpha,beta,mu)

  end function saiii3_efold_primitive


  
!  returns x at bfold=-efolds before the end of inflation
  function saiii3_x_trajectory(bfold,xend,alpha,beta,mu)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,beta,mu
    real(kp) :: saiii3_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: saiiiData
    
    saiiiData%real1 = alpha
    saiiiData%real2 = beta
    saiiiData%real3 = mu
    saiiiData%real4 = -bfold + saiii3_efold_primitive(xend,alpha,beta,mu)

    saiii3_x_trajectory = zbrent(find_saiii_x_trajectory,xend,saiiiXBig,tolFind,saiiiData)

  end function saiii3_x_trajectory

  
end module saiii3sr
