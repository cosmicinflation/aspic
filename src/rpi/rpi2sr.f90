!slow-roll functions for the R-R^p inflation potential
!
!V(phi) = M^4 * exp(-2y) * [ exp(y) - 1 ]^(2p/(2p-1))
!
!y = phi/Mp * sqrt(2/3)


module rpi2sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use rpicommon, only : rpi_norm_potential, rpi_norm_deriv_potential
  use rpicommon, only : rpi_norm_deriv_second_potential
  use rpicommon, only : rpi_epsilon_one, rpi_epsilon_two, rpi_epsilon_three
  use rpicommon, only : rpih_efold_primitive, rpih_x_trajectory, rpi_x_potmax
  use rpicommon, only : rpi_efold_primitive, find_rpi_x_trajectory
  implicit none

  private

  public rpi2_norm_potential, rpi2_epsilon_one, rpi2_epsilon_two, rpi2_epsilon_three
  public rpi2_efold_primitive, rpi2_x_trajectory
  public rpi2_norm_deriv_potential, rpi2_norm_deriv_second_potential

contains
 
  !returns V/M^4
  function rpi2_norm_potential(y,p)
    implicit none
    real(kp) :: rpi2_norm_potential
    real(kp), intent(in) :: y,p

    rpi2_norm_potential = rpi_norm_potential(y,p)

  end function rpi2_norm_potential



  !returns the first derivative of the potential with respect to y, divided by M^4
  function rpi2_norm_deriv_potential(y,p)
    implicit none
    real(kp) :: rpi2_norm_deriv_potential
    real(kp), intent(in) :: y,p

    rpi2_norm_deriv_potential = rpi_norm_deriv_potential(y,p)

  end function rpi2_norm_deriv_potential



  !returns the second derivative of the potential with respect to y, divided by M^4
  function rpi2_norm_deriv_second_potential(y,p)
    implicit none
    real(kp) :: rpi2_norm_deriv_second_potential
    real(kp), intent(in) :: y,p

    rpi2_norm_deriv_second_potential = rpi_norm_deriv_second_potential(y,p)

  end function rpi2_norm_deriv_second_potential



  !epsilon_one(y)
  function rpi2_epsilon_one(y,p)    
    implicit none
    real(kp) :: rpi2_epsilon_one
    real(kp), intent(in) :: y,p


    rpi2_epsilon_one = rpi_epsilon_one(y,p)


  end function rpi2_epsilon_one


  !epsilon_two(y)
  function rpi2_epsilon_two(y,p)    
    implicit none
    real(kp) :: rpi2_epsilon_two
    real(kp), intent(in) :: y,p

    rpi2_epsilon_two = rpi_epsilon_two(y,p)

  end function rpi2_epsilon_two


  !epsilon_three(y)
  function rpi2_epsilon_three(y,p)    
    implicit none
    real(kp) :: rpi2_epsilon_three
    real(kp), intent(in) :: y,p

    rpi2_epsilon_three = rpi_epsilon_three(y,p)

  end function rpi2_epsilon_three



  !this is integral[V(phi)/V'(phi) dphi]
  function rpi2_efold_primitive(y,p)
    implicit none
    real(kp), intent(in) :: y,p
    real(kp) :: rpi2_efold_primitive

    real(kp) :: yVmax

    if (p.eq.1._kp) then
       rpi2_efold_primitive = rpih_efold_primitive(y,p)
       return
    endif


    yVmax = rpi_x_potmax(p)

    if (y.lt.yVmax) stop 'rpi_efold_primitive: y < yVmax!'


    rpi2_efold_primitive = rpi_efold_primitive(y,p)
   
  end function rpi2_efold_primitive



  !returns y at bfold=-efolds before the end of inflation, ie N-Nend
  function rpi2_x_trajectory(bfold,yend,p)
    implicit none
    real(kp), intent(in) :: bfold, p, yend
    real(kp) :: rpi2_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rpi2Data
    real(kp) :: yVmax


    if (p.eq.1._kp) then !Higgs Inflation Model (HI)
       rpi2_x_trajectory = rpih_x_trajectory(bfold,yend,p)
       return
    endif

    yVmax = rpi_x_potmax(p)

    if (yend.lt.yVmax) stop 'rpi2_x_trajectory: yend < yVmax'
   
    mini = yVmax*(1._kp+epsilon(1._kp))
    maxi = yEnd

    rpi2Data%real1 = p
    rpi2Data%real2 = -bfold + rpi2_efold_primitive(yend,p)

    rpi2_x_trajectory = zbrent(find_rpi_x_trajectory,mini,maxi,tolFind,rpi2Data)

  end function rpi2_x_trajectory


end module rpi2sr
