!slow-roll functions for the running mass inflation 3 potential
!
!
!V(phi) = M^4 { 1 - c/2 [-1/2 + ln(x)](phi0*x)^2 }
!
!1: c<0 , x<1
!
!x = phi/phi0
!phi0 = phi0/Mp

module rmi3sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use rmicommon, only : rmi_norm_potential, rmi_norm_deriv_potential
  use rmicommon, only : rmi_norm_deriv_second_potential
  use rmicommon, only : rmi_epsilon_one, rmi_epsilon_two, rmi_epsilon_three
  use rmicommon, only : rmi_efold_primitive, find_rmi_x_trajectory


  implicit none

  private

  public rmi3_norm_potential, rmi3_norm_deriv_potential, rmi3_norm_deriv_second_potential
  public rmi3_epsilon_one, rmi3_epsilon_two, rmi3_epsilon_three
  public rmi3_efold_primitive, rmi3_x_trajectory
  public rmi3_numacc_xendmax

  
contains

!returns V/M**4
  function rmi3_norm_potential(x,c,phi0)    
    implicit none
    real(kp) :: rmi3_norm_potential
    real(kp), intent(in) :: x,c,phi0

    rmi3_norm_potential = rmi_norm_potential(x,c,phi0)

  end function rmi3_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function rmi3_norm_deriv_potential(x,c,phi0)
    implicit none
    real(kp) :: rmi3_norm_deriv_potential
    real(kp), intent(in) :: x,c,phi0
    
    rmi3_norm_deriv_potential = rmi_norm_deriv_potential(x,c,phi0)

  end function rmi3_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function rmi3_norm_deriv_second_potential(x,c,phi0)
    implicit none
    real(kp) :: rmi3_norm_deriv_second_potential
    real(kp), intent(in) :: x,c,phi0    

    rmi3_norm_deriv_second_potential &
         = rmi_norm_deriv_second_potential(x,c,phi0)
    

  end function rmi3_norm_deriv_second_potential



!epsilon_one(x)
  function rmi3_epsilon_one(x,c,phi0)    
    implicit none
    real(kp) :: rmi3_epsilon_one
    real(kp), intent(in) :: x,c,phi0

    rmi3_epsilon_one = rmi_epsilon_one(x,c,phi0)
    
  end function rmi3_epsilon_one


!epsilon_two(x)
  function rmi3_epsilon_two(x,c,phi0)    
    implicit none
    real(kp) :: rmi3_epsilon_two
    real(kp), intent(in) :: x,c,phi0

    rmi3_epsilon_two = rmi_epsilon_two(x,c,phi0)
    
  end function rmi3_epsilon_two


!epsilon_three(x)
  function rmi3_epsilon_three(x,c,phi0)    
    implicit none
    real(kp) :: rmi3_epsilon_three
    real(kp), intent(in) :: x,c,phi0
   
    rmi3_epsilon_three = rmi_epsilon_three(x,c,phi0)
    
  end function rmi3_epsilon_three


!this is integral[V(phi)/V'(phi) dphi]
  function rmi3_efold_primitive(x,c,phi0)
    implicit none
    real(kp), intent(in) :: x,c,phi0
    real(kp) :: rmi3_efold_primitive

    rmi3_efold_primitive = rmi_efold_primitive(x,c,phi0)
  
  end function rmi3_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function rmi3_x_trajectory(bfold,xend,c,phi0)
    implicit none
    real(kp), intent(in) :: bfold, xend, c, phi0
    real(kp) :: rmi3_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rmi3Data


    mini = epsilon(1._kp)
    maxi = xend*(1._kp-epsilon(1._kp))
  
    rmi3Data%real1 = c
    rmi3Data%real2 = phi0
    rmi3Data%real3 = -bfold + rmi3_efold_primitive(xend,c,phi0)
    
    rmi3_x_trajectory = zbrent(find_rmi_x_trajectory,mini,maxi,tolFind,rmi3Data)
       
  end function rmi3_x_trajectory


! Return an upper bound on xend for numerica computability
!  function rmi3_numacc_xendmax(c,phi0)
!    implicit none
!    real(kp), intent(in) :: c,phi0
!    real(kp) :: rmi3_numacc_xendmax

!Using an asymptotic expression for eps1 when x->1, and requiring
!eps1>epsilon(1._kp) for numerical convergence
!    rmi3_numacc_xendmax = 1._kp-sqrt(2._kp*epsilon(1._kp)* &
!         (1._kp+c*phi0**2/4._kp)**2/(c**2*phi0**2)) 


!  end function rmi3_numacc_xendmax
  
 function rmi3_numacc_xendmax(efold,c,phi0)
    implicit none
    real(kp), intent(in) :: efold,c,phi0
    real(kp) :: rmi3_numacc_xendmax,xMin, xMax, efoldMax
    real(kp), parameter :: tolFind=tolkp

 !Using an asymptotic expression for eps1 when x->1, and requiring
 !eps1>epsilon(1._kp) for numerical convergence

    xMin = sqrt(2._kp*epsilon(1._kp)/(c**2*phi0**2))

    xMax = 1._kp-sqrt(2._kp*epsilon(1._kp)* &
         (1._kp+c*phi0**2/4._kp)**2/(c**2*phi0**2)) 

    efoldMax = rmi3_efold_primitive(xMin,c,phi0) &
         - rmi3_efold_primitive(xMax,c,phi0)

    if (efold.gt.efoldMax) then
       write(*,*)'rmi3_xendmax: not enough efolds computable at current accuracy!'
       write(*,*)'efold requested=',efold,'   efold maxi=',efoldMax
       stop
    endif

    rmi3_numacc_xendmax = rmi3_x_trajectory(efold,xMax,c,phi0)
       

  end function rmi3_numacc_xendmax

end module rmi3sr
