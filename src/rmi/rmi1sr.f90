!slow-roll functions for the running mass inflation 1 potential
!
!
!V(phi) = M^4 { 1 - c/2 [-1/2 + ln(x)](phi0*x)^2 }
!
!1: c>0 , x<1
!
!x = phi/phi0
!phi0 = phi0/Mp

module rmi1sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use rmicommon, only : rmi_norm_potential, rmi_norm_deriv_potential,rmi_norm_deriv_second_potential
  use rmicommon, only : rmi_epsilon_one, rmi_epsilon_two, rmi_epsilon_three
  use rmicommon, only : rmi_efold_primitive, find_rmitraj


  implicit none

  private

  public rmi1_norm_potential, rmi1_norm_deriv_potential, rmi1_norm_deriv_second_potential
  public rmi1_epsilon_one, rmi1_epsilon_two, rmi1_epsilon_three
  public rmi1_efold_primitive, rmi1_x_trajectory, rmi1_xendmax
  
contains

!returns V/M**4
  function rmi1_norm_potential(x,c,phi0)    
    implicit none
    real(kp) :: rmi1_norm_potential
    real(kp), intent(in) :: x,c,phi0

    rmi1_norm_potential = rmi_norm_potential(x,c,phi0)

  end function rmi1_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function rmi1_norm_deriv_potential(x,c,phi0)
    implicit none
    real(kp) :: rmi1_norm_deriv_potential
    real(kp), intent(in) :: x,c,phi0
    
    rmi1_norm_deriv_potential = rmi_norm_deriv_potential(x,c,phi0)

  end function rmi1_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function rmi1_norm_deriv_second_potential(x,c,phi0)
    implicit none
    real(kp) :: rmi1_norm_deriv_second_potential
    real(kp), intent(in) :: x,c,phi0    

    rmi1_norm_deriv_second_potential &
         = rmi_norm_deriv_second_potential(x,c,phi0)
    

  end function rmi1_norm_deriv_second_potential



!epsilon_one(x)
  function rmi1_epsilon_one(x,c,phi0)    
    implicit none
    real(kp) :: rmi1_epsilon_one
    real(kp), intent(in) :: x,c,phi0

    rmi1_epsilon_one = rmi_epsilon_one(x,c,phi0)
    
  end function rmi1_epsilon_one


!epsilon_two(x)
  function rmi1_epsilon_two(x,c,phi0)    
    implicit none
    real(kp) :: rmi1_epsilon_two
    real(kp), intent(in) :: x,c,phi0

    rmi1_epsilon_two = rmi_epsilon_two(x,c,phi0)
    
  end function rmi1_epsilon_two


!epsilon_three(x)
  function rmi1_epsilon_three(x,c,phi0)    
    implicit none
    real(kp) :: rmi1_epsilon_three
    real(kp), intent(in) :: x,c,phi0
   
    rmi1_epsilon_three = rmi_epsilon_three(x,c,phi0)
    
  end function rmi1_epsilon_three


!this is integral[V(phi)/V'(phi) dphi]
  function rmi1_efold_primitive(x,c,phi0)
    implicit none
    real(kp), intent(in) :: x,c,phi0
    real(kp) :: rmi1_efold_primitive

    rmi1_efold_primitive = rmi_efold_primitive(x,c,phi0)
  
  end function rmi1_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function rmi1_x_trajectory(bfold,xend,c,phi0)
    implicit none
    real(kp), intent(in) :: bfold, xend, c, phi0
    real(kp) :: rmi1_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rmi1Data

    if (bfold .lt. 0._kp) then
    mini = xend*(1._kp+epsilon(1._kp))
    maxi = 1._kp*(1._kp-epsilon(1._kp))
    else
    mini = epsilon(1._kp)
    maxi = xend!*(1._kp-epsilon(1._kp))
    endif
  
    rmi1Data%real1 = c
    rmi1Data%real2 = phi0
    rmi1Data%real3 = -bfold + rmi1_efold_primitive(xend,c,phi0)
    
    rmi1_x_trajectory = zbrent(find_rmitraj,mini,maxi,tolFind,rmi1Data)
       
  end function rmi1_x_trajectory

!returns the maximal value for xend such that there are efold number
!of inflation from xtopNUM
  function rmi1_xendmax(efold,c,phi0)
    implicit none
    real(kp), intent(in) :: efold,c,phi0
    real(kp) :: rmi1_xendmax,xMin, xMax, efoldMax, eps
    real(kp), parameter :: tolFind=tolkp
   
    xMax =  1._kp-sqrt(2._kp*epsilon(1._kp)*(1._kp+c*phi0**2/4._kp)**2/(c**2*phi0**2)) !Using an asymptotic expression for eps1 when x->1, and requiring eps1>epsilon(1._kp) for numerical convergence
    xMin = sqrt(2._kp*epsilon(1._kp)/(c**2*phi0**2)) !Using an asymptotic expression for eps1 when x->0, and requiring eps1>epsilon(1._kp) for numerical convergence

    efoldMax = -rmi1_efold_primitive(xMin,c,phi0) &
         + rmi1_efold_primitive(xMax,c,phi0)

    if (efold.gt.efoldMax) then
       write(*,*)'rmi1_xendmax: not enough efolds!'
       write(*,*)'efold requested=',efold,'   efold maxi=',efoldMax
       stop
    endif

    rmi1_xendmax = rmi1_x_trajectory(efold,xMax,c,phi0)


  end function rmi1_xendmax
 

end module rmi1sr
