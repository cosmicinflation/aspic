!Common slow-roll functions for the running mass inflation 1,2,3,4 potential
!
!
!V(phi) = M**4 { 1 - c/2 [-1/2 + ln(x)](phi0*x)**2 }
!
!1: c>0 , x<1
!2: c>0 , x>1
!3: c<0 , x<1
!4: c<0 , x>1
!
!x = phi/phi0
!phi0 = phi0/Mp


module rmicommon
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : ei
  implicit none

  private

  public rmi_norm_potential
  public rmi_norm_deriv_potential, rmi_norm_deriv_second_potential
  public rmi_epsilon_one, rmi_epsilon_two, rmi_epsilon_three
  public rmi_efold_primitive, find_rmi_x_trajectory

contains


!returns V/M**4
  function rmi_norm_potential(x,c,phi0)
    implicit none
    real(kp) :: rmi_norm_potential
    real(kp), intent(in) :: x,c,phi0

    rmi_norm_potential = 1._kp+c/2._kp*(0.5_kp-log(x))*(phi0*x)**2

  end function rmi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function rmi_norm_deriv_potential(x,c,phi0)
    implicit none
    real(kp) :: rmi_norm_deriv_potential
    real(kp), intent(in) :: x,c,phi0

    rmi_norm_deriv_potential = -c*phi0**2*x*log(x)

  end function rmi_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function rmi_norm_deriv_second_potential(x,c,phi0)
    implicit none
    real(kp) :: rmi_norm_deriv_second_potential
    real(kp), intent(in) :: x,c,phi0

   
    rmi_norm_deriv_second_potential = -c*phi0**2*(1._kp+log(x))
    

  end function rmi_norm_deriv_second_potential



!epsilon_one(x)
  function rmi_epsilon_one(x,c,phi0)    
    implicit none
    real(kp) :: rmi_epsilon_one
    real(kp), intent(in) :: x,c,phi0


    rmi_epsilon_one = c**2/2._kp*(phi0*x)**2*log(x)**2/ &
         ((1._kp+c/2._kp*(0.5_kp-log(x))*(phi0*x)**2)**2)
    
  end function rmi_epsilon_one


!epsilon_two(x)
  function rmi_epsilon_two(x,c,phi0)    
    implicit none
    real(kp) :: rmi_epsilon_two
    real(kp), intent(in) :: x,c,phi0

    
    rmi_epsilon_two = (8._kp*c*(4._kp+c*(phi0*x)**2+ &
         log(x)*(4._kp-c*(phi0*x)**2+2._kp*c*(phi0*x)**2*log(x))))/ &
         (4._kp+c*(phi0*x)**2-2._kp*c*(phi0*x)**2*log(x))**2
   
  end function rmi_epsilon_two


!epsilon_three(x)
  function rmi_epsilon_three(x,c,phi0)    
    implicit none
    real(kp) :: rmi_epsilon_three
    real(kp), intent(in) :: x,c,phi0

    
    rmi_epsilon_three = (4._kp*c*log(x)*((4._kp+c*(phi0*x)**2)**2+ &
         8._kp*c*(phi0*x)**2*log(x)*(4._kp+c*(phi0*x)**2+ &
         log(x)*(6._kp-c*(phi0*x)**2+c*(phi0*x)**2*log(x)))))/((4._kp+ &
         c*(phi0*x)**2-2._kp*c*(phi0*x)**2*log(x))**2*(4._kp+c*(phi0*x)**2+ &
         log(x)*(4._kp-c*(phi0*x)**2+2._kp*c*(phi0*x)**2._kp*log(x))))
    
  end function rmi_epsilon_three


!this is integral[V(phi)/V'(phi) dphi]
  function rmi_efold_primitive(x,c,phi0)
    implicit none
    real(kp), intent(in) :: x,c,phi0
    real(kp) :: rmi_efold_primitive


    if (phi0 .eq. 0._kp) stop 'rmi_efold_primitive: phi/Mp=0!'

    rmi_efold_primitive = -log(abs(log(x)))/c+ &
         x**2*phi0**2/4._kp- &
         phi0**2*ei(2._kp*log(x))/4._kp

  end function rmi_efold_primitive



  function find_rmi_x_trajectory(x,rmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: rmiData
    real(kp) :: find_rmi_x_trajectory
    real(kp) :: c,phi0,NplusNuend

    c = rmiData%real1
    phi0 = rmiData%real2
    NplusNuend = rmiData%real3

    find_rmi_x_trajectory = rmi_efold_primitive(x,c,phi0) - NplusNuend
   
  end function find_rmi_x_trajectory


end module rmicommon
