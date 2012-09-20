!slow-roll functions for the R-R^p inflation potential
!
!V(phi) = M^4 * exp(-2y) * [ exp(y) - 1 ]^(2p/(2p-1))
!
!y = phi/Mp * sqrt(2/3)


module rpi1sr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  use rpicommon, only : rpi_norm_potential, rpi_norm_deriv_potential
  use rpicommon, only : rpi_norm_deriv_second_potential
  use rpicommon, only : rpi_epsilon_one, rpi_epsilon_two, rpi_epsilon_three
  use rpicommon, only : rpih_efold_primitive, rpih_x_trajectory, rpi_x_potmax
  implicit none

  private

  public rpi1_norm_potential, rpi1_epsilon_one, rpi1_epsilon_two, rpi1_epsilon_three
  public rpi1_x_endinf, rpi1_efold_primitive, rpi1_x_trajectory
  public rpi1_norm_deriv_potential, rpi1_norm_deriv_second_potential 

contains
  !returns V/M^4
  function rpi1_norm_potential(y,p)
    implicit none
    real(kp) :: rpi1_norm_potential
    real(kp), intent(in) :: y,p

    rpi1_norm_potential = rpi_norm_potential(y,p)

  end function rpi1_norm_potential



  !returns the first derivative of the potential with respect to y, divided by M^4
  function rpi1_norm_deriv_potential(y,p)
    implicit none
    real(kp) :: rpi1_norm_deriv_potential
    real(kp), intent(in) :: y,p

    rpi1_norm_deriv_potential = rpi_norm_deriv_potential(y,p)

  end function rpi1_norm_deriv_potential



  !returns the second derivative of the potential with respect to y, divided by M^4
  function rpi1_norm_deriv_second_potential(y,p)
    implicit none
    real(kp) :: rpi1_norm_deriv_second_potential
    real(kp), intent(in) :: y,p

    rpi1_norm_deriv_second_potential = rpi_norm_deriv_second_potential(y,p)

  end function rpi1_norm_deriv_second_potential



  !epsilon_one(y)
  function rpi1_epsilon_one(y,p)    
    implicit none
    real(kp) :: rpi1_epsilon_one
    real(kp), intent(in) :: y,p


    rpi1_epsilon_one = rpi_epsilon_one(y,p)


  end function rpi1_epsilon_one


  !epsilon_two(y)
  function rpi1_epsilon_two(y,p)    
    implicit none
    real(kp) :: rpi1_epsilon_two
    real(kp), intent(in) :: y,p

    rpi1_epsilon_two = rpi_epsilon_two(y,p)

  end function rpi1_epsilon_two


  !epsilon_three(y)
  function rpi1_epsilon_three(y,p)    
    implicit none
    real(kp) :: rpi1_epsilon_three
    real(kp), intent(in) :: y,p

    rpi1_epsilon_three = rpi_epsilon_three(y,p)

  end function rpi1_epsilon_three


 
  function rpi1_x_endinf(p)
    implicit none
    real(kp), intent(in) :: p
    real(kp) :: rpi1_x_endinf


    rpi1_x_endinf = log((1._kp+sqrt(3._kp)/2._kp)/ &
         ((p-1._kp)/(2._kp*p-1._kp)+sqrt(3._kp)/2._kp))
  
  end function rpi1_x_endinf


 
  !this is integral[V(phi)/V'(phi) dphi]
  function rpi1_efold_primitive(y,p)
    implicit none
    real(kp), intent(in) :: y,p
    real(kp) :: rpi1_efold_primitive

    real(kp) :: yVmax

    if (p.eq.1._kp) then
       rpi1_efold_primitive = rpih_efold_primitive(y,p)
       return
    endif


    yVmax = rpi_x_potmax(p)

    if (y.gt.yVmax) stop 'rpi1_efold_primitive: y > yVmax'

    rpi1_efold_primitive = 3._kp/4._kp*(-p/(p-1._kp)*log(-exp(y)- &
         (1._kp-2._kp*p)/(p-1._kp))-y)

  end function rpi1_efold_primitive


  !returns y at bfold=-efolds before the end of inflation, ie N-Nend
  function rpi1_x_trajectory(bfold,yend,p)
    implicit none
    real(kp), intent(in) :: bfold, p, yend
    real(kp) :: rpi1_x_trajectory
    real(kp) :: yVmax


    if (p.eq.1._kp) then !Higgs Inflation Model (HI)
       rpi1_x_trajectory = rpih_x_trajectory(bfold,yend,p)
       return
    endif

    yVmax = rpi_x_potmax(p)

    if (yend.gt.yVmax) stop 'rpi1_x_trajectory: yend > yVmax'

    rpi1_x_trajectory = log(-p/(p-1._kp)*(lambert(-((1._kp-2._kp*p)/p+(p-1._kp)*exp(yend)/p)* &
         exp(-(1._kp-2._kp)/p-(p-1._kp)/p*exp(yend)+3._kp/4._kp*(p-1._kp)/p*bfold),0) &
         -(2._kp*p-1._kp)/p))


  end function rpi1_x_trajectory


end module rpi1sr
