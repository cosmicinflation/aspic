!slow-roll functions for the R-R^p inflation potential
!
!V(phi) = M^4 * exp(-2y) * [ exp(y) - 1 ]^(2p/(2p-1))
!
!y = phi/Mp * sqrt(2/3)


module rpisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  implicit none

  private

  public  rpi_norm_potential, rpi_epsilon_one, rpi_epsilon_two, rpi_epsilon_three
  public  rpi_y_endinf, rpi_efold_primitive, rpi_y_trajectory
  public  rpi_norm_deriv_potential, rpi_norm_deriv_second_potential
 
contains
!returns V/M^4
  function rpi_norm_potential(y,p)
    implicit none
    real(kp) :: rpi_norm_potential
    real(kp), intent(in) :: y,p

    rpi_norm_potential = exp(-2._kp*y)*(exp(y)-1._kp)**(2._kp*p/(2._kp*p-1._kp))

  end function rpi_norm_potential



!returns the first derivative of the potential with respect to x=phi/Mp, divided by M^4
  function rpi_norm_deriv_potential(y,p)
    implicit none
    real(kp) :: rpi_norm_deriv_potential
    real(kp), intent(in) :: y,p

   rpi_norm_deriv_potential = (2._kp*sqrt(2._kp/3._kp)*exp(-2._kp*y)* & 
                               (-1._kp+exp(y))**(1._kp/(-1._kp+2._kp*p)) * &
                               (1._kp+exp(y)*(-1._kp+p)-2._kp*p))/(1._kp-2._kp*p)

  end function rpi_norm_deriv_potential



!returns the second derivative of the potential with respect to x=phi/Mp, divided by M^4
  function rpi_norm_deriv_second_potential(y,p)
    implicit none
    real(kp) :: rpi_norm_deriv_second_potential
    real(kp), intent(in) :: y,p

    rpi_norm_deriv_second_potential = 2._kp/((1._kp-2._kp*p)**2)*sqrt(2._kp/3._kp)* &
                                      exp(-y)*(-1._kp+exp(y))**(-1._kp+1._kp/(-1._kp+2._kp*p))* &
                                      (-4._kp+(13._kp-10._kp*p)*p+2._kp*(2._kp+p*(-6._kp+5._kp*p))* &
                                      cosh(y)+2._kp*(2._kp-3._kp*p)*p*sinh(y))

  end function rpi_norm_deriv_second_potential



!epsilon_one(y)
  function rpi_epsilon_one(y,p)    
    implicit none
    real(kp) :: rpi_epsilon_one
    real(kp), intent(in) :: y,p

  
    rpi_epsilon_one = (4._kp*(1._kp+exp(y)*(-1._kp+p)-2._kp*p)**2)/ &
                      (3._kp*(-1._kp+exp(y))**2*(1._kp-2._kp*p)**2)

    
  end function rpi_epsilon_one


!epsilon_two(y)
  function rpi_epsilon_two(y,p)    
    implicit none
    real(kp) :: rpi_epsilon_two
    real(kp), intent(in) :: y,p
    
    rpi_epsilon_two = (2._kp*p*(1._kp/sinh(y/2._kp))**2)/(-3._kp+6._kp*p)
    
  end function rpi_epsilon_two


!epsilon_three(y)
  function rpi_epsilon_three(y,p)    
    implicit none
    real(kp) :: rpi_epsilon_three
    real(kp), intent(in) :: y,p
    
    rpi_epsilon_three = -((4._kp*(1._kp+exp(y))*(1._kp+exp(y)*(-1._kp+p)-2._kp*p))/ &
                         (3._kp*(-1._kp+exp(y))**2*(-1._kp+2._kp*p)))
    
  end function rpi_epsilon_three


!returns y at the end of inflation defined as epsilon1=1
  function rpi_y_endinf(p)
    implicit none
    real(kp), intent(in) :: p
    real(kp) :: rpi_y_endinf
    

    rpi_y_endinf = log((1._kp+sqrt(3._kp)/2._kp)/ &
                   ((p-1._kp)/(2._kp*p-1._kp)+sqrt(3._kp)/2._kp))


  end function rpi_y_endinf




!this is integral[V(phi)/V'(phi) dphi]
  function rpi_efold_primitive(y,p)
    implicit none
    real(kp), intent(in) :: y,p
    real(kp) :: rpi_efold_primitive

    if (p.eq.1._kp) then

	rpi_efold_primitive = 3._kp/4._kp*(exp(y)-y) !Higgs Inflation Model (HI)
  
    else

        !If inflation proceeds from the left to the right
!        rpi_efold_primitive = 3._kp/4._kp*(-p/(p-1._kp)*log(exp(y)+ &
!                             (1._kp-2._kp*p)/(p-1._kp))-y)

        !If inflation proceeds from the right to the left
        rpi_efold_primitive = 3._kp/4._kp*(-p/(p-1._kp)*log(-exp(y)- &
                              (1._kp-2._kp*p)/(p-1._kp))-y)

    endif

  end function rpi_efold_primitive


!returns y at bfold=-efolds before the end of inflation, ie N-Nend
  function rpi_y_trajectory(bfold,yend,p)
    implicit none
    real(kp), intent(in) :: bfold, p, yend
    real(kp) :: rpi_y_trajectory


    if (p.eq.1._kp) then !Higgs Inflation Model (HI)

	rpi_y_trajectory = (4._kp/3._kp*bfold+yend-exp(yend)- &
                           lambert(-exp(4._kp/3._kp*bfold+yend-exp(yend)),-1))

    else

    ! If Inflation proceeds from the left to the right 
!    rpi_y_trajectory = log(p/(p-1._kp)*(lambert(((1._kp-2._kp*p)/p+(p-1._kp)*exp(yend)/p)* &
!                       exp((1._kp-2._kp)/p(p-1._kp)/p*exp(yend)+3._kp/4._kp*(p-1._kp)/p*bfold),0) &
!                       +(2._kp-1._kp)/p))


    ! If Inflation proceeds from the right to the left 
    rpi_y_trajectory = log(-p/(p-1._kp)*(lambert(-((1._kp-2._kp*p)/p+(p-1._kp)*exp(yend)/p)* &
                       exp(-(1._kp-2._kp)/p-(p-1._kp)/p*exp(yend)+3._kp/4._kp*(p-1._kp)/p*bfold),0) &
                       -(2._kp-1._kp)/p))

    endif
       
  end function rpi_y_trajectory


end module rpisr
