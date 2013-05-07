!slow-roll functions for the open string tachyonic inflation potential
!
!V(phi) = - M**4 * x**2 Ln[x**2]
!
!x=phi/phi0
!phi0=phi0/Mp


module ostisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : ei
  use inftools, only : zbrent
  implicit none

  private

  public  osti_norm_potential, osti_epsilon_one, osti_epsilon_two, osti_epsilon_three
  public  osti_x_endinf, osti_efold_primitive, osti_x_trajectory
  public  osti_norm_deriv_potential, osti_norm_deriv_second_potential
 
contains
!returns V/M**4
  function osti_norm_potential(x,phi0)
    implicit none
    real(kp) :: osti_norm_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: phi0

    osti_norm_potential = -x**2*log(x**2)

  end function osti_norm_potential



!returns the first derivative of the potential with respect to x=phi/phi0, divided by M**4
  function osti_norm_deriv_potential(x,phi0)
    implicit none
    real(kp) :: osti_norm_deriv_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: phi0

   osti_norm_deriv_potential = -2._kp*x*(1._kp+log(x**2))

  end function osti_norm_deriv_potential



!returns the second derivative of the potential with respect to x=phi/phi0, divided by M**4
  function osti_norm_deriv_second_potential(x,phi0)
    implicit none
    real(kp) :: osti_norm_deriv_second_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: phi0

    osti_norm_deriv_second_potential = -2._kp*(3._kp+log(x**2))

  end function osti_norm_deriv_second_potential



!epsilon_one(x)
  function osti_epsilon_one(x,phi0)    
    implicit none
    real(kp) :: osti_epsilon_one
    real(kp), intent(in) :: x,phi0

  
    osti_epsilon_one = (2._kp*(1._kp+log(x**2))**2)/(phi0**2*x**2*log(x**2)**2)

    
  end function osti_epsilon_one


!epsilon_two(x)
  function osti_epsilon_two(x,phi0)    
    implicit none
    real(kp) :: osti_epsilon_two
    real(kp), intent(in) :: x,phi0
    
    osti_epsilon_two = (4._kp*(2._kp+log(x**2)+log(x**2)**2))/ &
                       (phi0**2*x**2*log(x**2)**2)
 
  end function osti_epsilon_two


!epsilon_three(x)
  function osti_epsilon_three(x,phi0)    
    implicit none
    real(kp) :: osti_epsilon_three
    real(kp), intent(in) :: x,phi0
    
    osti_epsilon_three = (4._kp*(1._kp+log(x**2))* &
                         (4._kp+log(x**2)*(3._kp+log(x**2)+log(x**2)**2)))/ &
                         (phi0**2*x**2*log(x**2)**2*(2._kp+log(x**2)+ &
                         log(x**2)**2))
   
  end function osti_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function osti_x_endinf(phi0)
    implicit none
    real(kp), intent(in) :: phi0
    real(kp) :: osti_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ostiData

    mini=epsilon(1._kp)
    maxi=exp(-0.5_kp)*(1._kp-epsilon(1._kp)) !Potential Maxiphi0m, where eps1 vanishes

    ostiData%real1=phi0

    osti_x_endinf = zbrent(find_osti_x_endinf,mini,maxi,tolFind,ostiData)


  end function osti_x_endinf


  function find_osti_x_endinf(x,ostiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ostiData
    real(kp) :: find_osti_x_endinf
    real(kp) :: phi0

    phi0 = ostiData%real1
    
    find_osti_x_endinf = osti_epsilon_one(x,phi0) - 1._kp
  
  end function find_osti_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function osti_efold_primitive(x,phi0)
    implicit none
    real(kp), intent(in) :: x,phi0
    real(kp) :: osti_efold_primitive


    osti_efold_primitive = 0.25_kp*phi0**2*(x**2-ei(1._kp+log(x**2))/exp(1._kp))


  end function osti_efold_primitive

!returns x=phi/mi at bfold=-efolds before the end of inflation, ie N-Nend
  function osti_x_trajectory(bfold,xend,phi0)
    implicit none
    real(kp), intent(in) :: bfold, xend, phi0
    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: osti_x_trajectory,mini,maxi
    type(transfert) :: ostiData

    ostiData%real1 = bfold
    ostiData%real2 = xend
    ostiData%real3 = phi0

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = exp(-0.5_kp)*(1._kp-epsilon(1._kp))
    
    osti_x_trajectory =zbrent(find_osti_x_trajectory,mini,maxi,tolzbrent,ostiData)
       
  end function osti_x_trajectory

  function find_osti_x_trajectory(x,ostiData)   
    implicit none
    real(kp) :: find_osti_x_trajectory
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ostiData

    real(kp) :: bfold,xend,phi0

    bfold = ostiData%real1
    xend = ostiData%real2
    phi0 = ostiData%real3

    find_osti_x_trajectory = osti_efold_primitive(x,phi0) - osti_efold_primitive(xend,phi0) + bfold
  
  end function find_osti_x_trajectory


end module ostisr
