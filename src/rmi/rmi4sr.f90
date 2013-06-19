!slow-roll functions for the running mass inflation 1 potential
!
!
!V(phi) = M^4 { 1 - c/2 [-1/2 + ln(x)](phi0*x)^2 }
!
!1: c<0 , x>1
!
!x = phi/phi0
!phi0 = phi0/Mp

module rmi4sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use rmicommon, only : rmi_norm_potential, rmi_norm_deriv_potential
  use rmicommon, only : rmi_norm_deriv_second_potential
  use rmicommon, only : rmi_epsilon_one, rmi_epsilon_two, rmi_epsilon_three
  use rmicommon, only : rmi_efold_primitive, find_rmi_x_trajectory


  implicit none

  private

  public rmi4_norm_potential, rmi4_norm_deriv_potential, rmi4_norm_deriv_second_potential
  public rmi4_epsilon_one, rmi4_epsilon_two, rmi4_epsilon_three
  public rmi4_efold_primitive, rmi4_x_trajectory, rmi4_numacc_xendmin

  
contains

!returns V/M**4
  function rmi4_norm_potential(x,c,phi0)    
    implicit none
    real(kp) :: rmi4_norm_potential
    real(kp), intent(in) :: x,c,phi0

    rmi4_norm_potential = rmi_norm_potential(x,c,phi0)

  end function rmi4_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function rmi4_norm_deriv_potential(x,c,phi0)
    implicit none
    real(kp) :: rmi4_norm_deriv_potential
    real(kp), intent(in) :: x,c,phi0
    
    rmi4_norm_deriv_potential = rmi_norm_deriv_potential(x,c,phi0)

  end function rmi4_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function rmi4_norm_deriv_second_potential(x,c,phi0)
    implicit none
    real(kp) :: rmi4_norm_deriv_second_potential
    real(kp), intent(in) :: x,c,phi0    

    rmi4_norm_deriv_second_potential &
         = rmi_norm_deriv_second_potential(x,c,phi0)
    

  end function rmi4_norm_deriv_second_potential



!epsilon_one(x)
  function rmi4_epsilon_one(x,c,phi0)    
    implicit none
    real(kp) :: rmi4_epsilon_one
    real(kp), intent(in) :: x,c,phi0

    rmi4_epsilon_one = rmi_epsilon_one(x,c,phi0)
    
  end function rmi4_epsilon_one


!epsilon_two(x)
  function rmi4_epsilon_two(x,c,phi0)    
    implicit none
    real(kp) :: rmi4_epsilon_two
    real(kp), intent(in) :: x,c,phi0

    rmi4_epsilon_two = rmi_epsilon_two(x,c,phi0)
    
  end function rmi4_epsilon_two


!epsilon_three(x)
  function rmi4_epsilon_three(x,c,phi0)    
    implicit none
    real(kp) :: rmi4_epsilon_three
    real(kp), intent(in) :: x,c,phi0
   
    rmi4_epsilon_three = rmi_epsilon_three(x,c,phi0)
    
  end function rmi4_epsilon_three


!this is integral[V(phi)/V'(phi) dphi]
  function rmi4_efold_primitive(x,c,phi0)
    implicit none
    real(kp), intent(in) :: x,c,phi0
    real(kp) :: rmi4_efold_primitive

    rmi4_efold_primitive = rmi_efold_primitive(x,c,phi0)
  
  end function rmi4_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function rmi4_x_trajectory(bfold,xend,c,phi0)
    implicit none
    real(kp), intent(in) :: bfold, xend, c, phi0
    real(kp) :: rmi4_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rmi4Data


    mini = xend*(1._kp+epsilon(1._kp))
    maxi = 1._kp/epsilon(1._kp)
  
    rmi4Data%real1 = c
    rmi4Data%real2 = phi0
    rmi4Data%real3 = -bfold + rmi4_efold_primitive(xend,c,phi0)
    
    rmi4_x_trajectory = zbrent(find_rmi_x_trajectory,mini,maxi,tolFind,rmi4Data)
       
  end function rmi4_x_trajectory

! Return an upper bound on xend for numerica computability
  function rmi4_numacc_xendmin(c,phi0)
    implicit none
    real(kp), intent(in) :: c,phi0
    real(kp) :: rmi4_numacc_xendmin

 	rmi4_numacc_xendmin = 1._kp+sqrt(2._kp*epsilon(1._kp)* &
                              (1._kp+c*phi0**2/4._kp)**2/(c**2*phi0**2)) 
    !Using an asymptotic expression for eps1 when x->1, and requiring eps1>epsilon(1._kp) for numerical convergence


  end function rmi4_numacc_xendmin

end module rmi4sr
