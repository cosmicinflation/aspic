!slow-roll functions for the running mass inflation 2 potential
!
!
!V(phi) = M^4 { 1 - c/2 [-1/2 + ln(x)](phi0*x)^2 }
!
!1: c>0 , x>1
!
!x = phi/phi0
!phi0 = phi0/Mp

module rmi2sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use rmicommon, only : rmi_norm_potential, rmi_norm_deriv_potential
  use rmicommon, only : rmi_norm_deriv_second_potential
  use rmicommon, only : rmi_epsilon_one, rmi_epsilon_two, rmi_epsilon_three
  use rmicommon, only : rmi_efold_primitive, find_rmi_x_trajectory


  implicit none

  private

  public rmi2_norm_potential, rmi2_norm_deriv_potential, rmi2_norm_deriv_second_potential
  public rmi2_epsilon_one, rmi2_epsilon_two, rmi2_epsilon_three
  public rmi2_efold_primitive, rmi2_x_trajectory, rmi2_numacc_xendmin

contains

!returns V/M**4
  function rmi2_norm_potential(x,c,phi0)    
    implicit none
    real(kp) :: rmi2_norm_potential
    real(kp), intent(in) :: x,c,phi0

    rmi2_norm_potential = rmi_norm_potential(x,c,phi0)

  end function rmi2_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function rmi2_norm_deriv_potential(x,c,phi0)
    implicit none
    real(kp) :: rmi2_norm_deriv_potential
    real(kp), intent(in) :: x,c,phi0
    
    rmi2_norm_deriv_potential = rmi_norm_deriv_potential(x,c,phi0)

  end function rmi2_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function rmi2_norm_deriv_second_potential(x,c,phi0)
    implicit none
    real(kp) :: rmi2_norm_deriv_second_potential
    real(kp), intent(in) :: x,c,phi0    

    rmi2_norm_deriv_second_potential &
         = rmi_norm_deriv_second_potential(x,c,phi0)
    

  end function rmi2_norm_deriv_second_potential



!epsilon_one(x)
  function rmi2_epsilon_one(x,c,phi0)    
    implicit none
    real(kp) :: rmi2_epsilon_one
    real(kp), intent(in) :: x,c,phi0

    rmi2_epsilon_one = rmi_epsilon_one(x,c,phi0)
    
  end function rmi2_epsilon_one


!epsilon_two(x)
  function rmi2_epsilon_two(x,c,phi0)    
    implicit none
    real(kp) :: rmi2_epsilon_two
    real(kp), intent(in) :: x,c,phi0

    rmi2_epsilon_two = rmi_epsilon_two(x,c,phi0)
    
  end function rmi2_epsilon_two


!epsilon_three(x)
  function rmi2_epsilon_three(x,c,phi0)    
    implicit none
    real(kp) :: rmi2_epsilon_three
    real(kp), intent(in) :: x,c,phi0
   
    rmi2_epsilon_three = rmi_epsilon_three(x,c,phi0)
    
  end function rmi2_epsilon_three


!this is integral[V(phi)/V'(phi) dphi]
  function rmi2_efold_primitive(x,c,phi0)
    implicit none
    real(kp), intent(in) :: x,c,phi0
    real(kp) :: rmi2_efold_primitive

    rmi2_efold_primitive = rmi_efold_primitive(x,c,phi0)
  
  end function rmi2_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function rmi2_x_trajectory(bfold,xend,c,phi0)
    implicit none
    real(kp), intent(in) :: bfold, xend, c, phi0
    real(kp) :: rmi2_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rmi2Data
    
    if (bfold .lt. 0._kp) then
       mini = 1._kp*(1._kp+epsilon(1._kp))
       maxi = xend*(1._kp-epsilon(1._kp))
    else
       mini = xend
       maxi = 2._kp
    endif
  
    rmi2Data%real1 = c
    rmi2Data%real2 = phi0
    rmi2Data%real3 = -bfold + rmi2_efold_primitive(xend,c,phi0)
    
    rmi2_x_trajectory = zbrent(find_rmi_x_trajectory,mini,maxi,tolFind,rmi2Data)
       
  end function rmi2_x_trajectory

!returns the minimal value for xend such that there are efold number
!of inflation from xtopNUM
  function rmi2_numacc_xendmin(efold,c,phi0)
    implicit none
    real(kp), intent(in) :: efold,c,phi0
    real(kp) :: rmi2_numacc_xendmin,xMin, xMax, efoldMax, eps
    real(kp), parameter :: tolFind=tolkp

 !Using an asymptotic expression for eps1 when x->1, and requiring
 !eps1>epsilon(1._kp) for numerical convergence
    xMin =  1._kp+sqrt(2._kp*epsilon(1._kp)*(1._kp+c*phi0**2/4._kp)**2/(c**2*phi0**2))
    xMax = 2._kp

    efoldMax = rmi2_efold_primitive(xMin,c,phi0) &
         - rmi2_efold_primitive(xMax,c,phi0)

    if (efold.gt.efoldMax) then
       write(*,*)'rmi2_xendmax: not enough efolds computable at current accuracy!'
       write(*,*)'efold requested=',efold,'   efold maxi=',efoldMax
       stop
    endif

    rmi2_numacc_xendmin = rmi2_x_trajectory(efold,xMin,c,phi0)
       

  end function rmi2_numacc_xendmin
 

end module rmi2sr
