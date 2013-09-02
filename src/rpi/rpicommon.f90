!slow-roll functions for the R-R^p inflation potential
!
!V(phi) = M^4 * exp(-2y) * [ exp(y) - 1 ]^(2p/(2p-1))
!
!y = phi/Mp * sqrt(2/3)


module rpicommon
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  use inftools, only : zbrent
  implicit none

  private

  public rpi_norm_potential, rpi_epsilon_one, rpi_epsilon_two, rpi_epsilon_three
  public rpi_norm_deriv_potential, rpi_norm_deriv_second_potential
  public rpih_x_trajectory, rpih_efold_primitive, rpi_x_potmax
  public rpi_x_epsoneunity, rpi_efold_primitive, find_rpi_x_trajectory

!because we use exp(x) everywhere
  real(kp), parameter :: xBig = log(epsilon(1._kp)*huge(1._kp))

  public xBig

contains
  !returns V/M^4
  function rpi_norm_potential(y,p)
    implicit none
    real(kp) :: rpi_norm_potential
    real(kp), intent(in) :: y,p

!    rpi_norm_potential = exp(-2._kp*y)*(exp(y)-1._kp)**(2._kp*p/(2._kp*p-1._kp))

    rpi_norm_potential = ((exp(-(p-1._kp)*y/p) &
         - exp(-(2._kp*p-1._kp)*y/p))**2)**(p/(2._kp*p-1))

  end function rpi_norm_potential



!returns the first derivative of the potential with respect to y, divided by M^4
  function rpi_norm_deriv_potential(y,p)
    implicit none
    real(kp) :: rpi_norm_deriv_potential
    real(kp), intent(in) :: y,p

!    rpi_norm_deriv_potential = (2._kp*exp(-2._kp*y)* & 
!         (-1._kp+exp(y))**(1._kp/(-1._kp+2._kp*p)) * &
!         (1._kp+exp(y)*(-1._kp+p)-2._kp*p))/(1._kp-2._kp*p)

    rpi_norm_deriv_potential = ( 2._kp &
         * (exp(-(2._kp*p-2._kp)*y) - exp(-(2._kp*p-1._kp)*y))**(1._kp/(-1._kp+2._kp*p)) &
         *  ( (1._kp-2._kp*p)*exp(-y) + p -1._kp ) )/(1._kp-2._kp*p)

  end function rpi_norm_deriv_potential



  !returns the second derivative of the potential with respect to y, divided by M^4
  function rpi_norm_deriv_second_potential(y,p)
    implicit none
    real(kp) :: rpi_norm_deriv_second_potential
    real(kp), intent(in) :: y,p

!    rpi_norm_deriv_second_potential = 2._kp/((1._kp-2._kp*p)**2)* &
!         exp(-y)*(-1._kp+exp(y))**(-1._kp+1._kp/(-1._kp+2._kp*p))* &
!         (-4._kp+(13._kp-10._kp*p)*p+2._kp*(2._kp+p*(-6._kp+5._kp*p))* &
!         cosh(y)+2._kp*(2._kp-3._kp*p)*p*sinh(y))
  
    rpi_norm_deriv_second_potential = 2._kp/((1._kp-2._kp*p)**2)* &
         exp(-y*(4._kp*p-3._kp)/(2._kp*p-1._kp)) &
         *(1._kp- exp(-y))**(-1._kp+1._kp/(-1._kp+2._kp*p))* &
         (-4._kp+(13._kp-10._kp*p)*p+2._kp*(2._kp+p*(-6._kp+5._kp*p))* &
         cosh(y)+2._kp*(2._kp-3._kp*p)*p*sinh(y))

  end function rpi_norm_deriv_second_potential



  !epsilon_one(y)
  function rpi_epsilon_one(y,p)    
    implicit none
    real(kp) :: rpi_epsilon_one
    real(kp), intent(in) :: y,p

! exp(y) is problematic at infinity ;)
!    rpi_epsilon_one = (4._kp*(1._kp+exp(y)*(-1._kp+p)-2._kp*p)**2)/ &
!         (3._kp*(-1._kp+exp(y))**2*(1._kp-2._kp*p)**2)

    rpi_epsilon_one = (4._kp*((1._kp - 2._kp*p)*exp(-y) +(-1._kp+p))**2)/ &
         (3._kp*(-exp(-y)+1._kp)**2*(1._kp-2._kp*p)**2)
    

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

!    rpi_epsilon_three = -((4._kp*(1._kp+exp(y))*(1._kp+exp(y)*(-1._kp+p)-2._kp*p))/ &
!         (3._kp*(-1._kp+exp(y))**2*(-1._kp+2._kp*p)))

    rpi_epsilon_three = -((4._kp*(1._kp+exp(-y))*((1._kp-2._kp*p)*exp(-y) &
         +(-1._kp+p))) / &
         (3._kp*(1._kp-exp(-y))**2*(-1._kp+2._kp*p)))
 

  end function rpi_epsilon_three


!returns the solution of eps1=1, which is xend for rpi1 and rpi3
  function rpi_x_epsoneunity(p)
    implicit none
    real(kp), intent(in) :: p
    real(kp) :: rpi_x_epsoneunity

    if (p.eq.1._kp) then

       rpi_x_epsoneunity = sqrt(3._kp/2._kp)*log(1._kp+2._kp/sqrt(3._kp))

    else

       rpi_x_epsoneunity = log((1._kp+sqrt(3._kp)/2._kp)/ &
            ((p-1._kp)/(2._kp*p-1._kp)+sqrt(3._kp)/2._kp))

    endif

  end function rpi_x_epsoneunity


!field value at which the potential is maximal
  function rpi_x_potmax(p)
    implicit none
    real(kp) , intent(in) :: p
    real(kp) :: rpi_x_potmax

    if (p.eq.1._kp) stop 'rpi_x_potmax: no maximum, p=1 is singular!'

    rpi_x_potmax = log((2._kp*p-1._kp)/(p-1._kp))

  end function rpi_x_potmax

!this is integral[V(phi)/V'(phi) dphi]
  function rpi_efold_primitive(y,p)
    implicit none
    real(kp), intent(in) :: y,p
    real(kp) :: rpi_efold_primitive

    if (p.eq.1._kp) then
       rpi_efold_primitive = rpih_efold_primitive(y,p)
       return
    endif

!    rpi_efold_primitive = 3._kp/4._kp*(-p/(p-1._kp)*log(abs(-exp(y)- &
!         (1._kp-2._kp*p)/(p-1._kp)))-y)

    rpi_efold_primitive = 3._kp/4._kp*( -p/(p-1._kp)*log(abs(1._kp + &
         exp(-y)*(1._kp-2._kp*p)/(p-1._kp))) + (1._kp-2._kp*p)/(p-1._kp)*y )

  end function rpi_efold_primitive

  function find_rpi_x_trajectory(y,rpiData)    
    implicit none
    real(kp), intent(in) :: y   
    type(transfert), optional, intent(inout) :: rpiData
    real(kp) :: find_rpi_x_trajectory
    real(kp) :: p,NplusNuend

    p = rpiData%real1
    NplusNuend = rpiData%real2

    find_rpi_x_trajectory = rpi_efold_primitive(y,p) - NplusNuend
   
  end function find_rpi_x_trajectory


!Higgs Inflation Model (HI)
  function rpih_efold_primitive(y,p)
    implicit none
    real(kp), intent(in) :: y,p
    real(kp) :: rpih_efold_primitive

    if (p.ne.1._kp) stop 'rpih_efold_primitive: p is not unity!'
   

    rpih_efold_primitive = 3._kp/4._kp*(exp(y)-y) 

  end function rpih_efold_primitive


  
  function rpih_x_trajectory(bfold,yend,p)
    implicit none
    real(kp), intent(in) :: bfold, p, yend
    real(kp) :: rpih_x_trajectory

    if (p.ne.1._kp) stop 'rpih_x_trajectory: p is not unity!'

    rpih_x_trajectory = (4._kp/3._kp*bfold+yend-exp(yend)- &
         lambert(-exp(4._kp/3._kp*bfold+yend-exp(yend)),-1))
      
  end function rpih_x_trajectory


end module rpicommon
