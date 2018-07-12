!slow-roll functions for the hybrid natural inflation potential
!
!V(phi) = M**4 [ 1 + alpha cos(x) )
!
!x = phi/phi0

module hnisr
  use infprec, only : kp, tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public hni_norm_potential, hni_norm_deriv_potential, hni_norm_deriv_second_potential
  public hni_epsilon_one, hni_epsilon_two, hni_epsilon_three
  public hni_efold_primitive, hni_x_trajectory, hni_x_endinf
  public hni_alphamin, hni_phizeromax


contains


  function hni_check_params(alpha,phi0)
    implicit none
    logical :: hni_check_params
    real(kp), intent(in) :: alpha,phi0

    hni_check_params = (alpha .ge. hni_alphamin(phi0))

  end function hni_check_params


!returns V/M**4
  function hni_norm_potential(x,alpha,phi0)
    implicit none
    real(kp) :: hni_norm_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: phi0

    hni_norm_potential = 1._kp+alpha*cos(x)

  end function hni_norm_potential


!returns the first derivative of the potential with respect to x=phi/phi0, divided by M**4
  function hni_norm_deriv_potential(x,alpha,phi0)
    implicit none
    real(kp) :: hni_norm_deriv_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: phi0

    hni_norm_deriv_potential = -alpha*sin(x)

  end function hni_norm_deriv_potential


!returns the second derivative of the potential with respect to x=phi/phi0, divided by M**4
  function hni_norm_deriv_second_potential(x,alpha,phi0)
    implicit none
    real(kp) :: hni_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: phi0

    hni_norm_deriv_second_potential = -alpha*cos(x)

  end function hni_norm_deriv_second_potential


!epsilon1(x)
  function hni_epsilon_one(x,alpha,phi0)    
    implicit none
    real(kp) :: hni_epsilon_one
    real(kp), intent(in) :: x,alpha,phi0

    hni_epsilon_one = ((alpha**2*sin(x)**2)/(2._kp*(1._kp+alpha*cos(x))**2))/phi0**2

  end function hni_epsilon_one


!epsilon2(x)
  function hni_epsilon_two(x,alpha,phi0)    
    implicit none
    real(kp) :: hni_epsilon_two
    real(kp), intent(in) :: x,alpha,phi0

    hni_epsilon_two = 2._kp/phi0**2*(alpha*(alpha+cos(x)))/(1._kp+alpha*cos(x))**2

  end function hni_epsilon_two

!epsilon3(x)
  function hni_epsilon_three(x,alpha,phi0)    
    implicit none
    real(kp) :: hni_epsilon_three
    real(kp), intent(in) :: x,alpha,phi0

    hni_epsilon_three = ((alpha*(-1._kp+2._kp*alpha**2+alpha*cos(x))* &
                        sin(x)**2)/((alpha+cos(x))*(1._kp+alpha*cos(x))**2))/phi0**2

  end function hni_epsilon_three

!returns x at the end of inflation defined as epsilon1=1
  function hni_x_endinf(alpha,phi0)
    implicit none
    real(kp) :: hni_x_endinf
    real(kp), intent(in) :: alpha,phi0

    hni_x_endinf =acos((sqrt(alpha**2/(4._kp*phi0**2)+alpha**2/2._kp-0.5_kp)/phi0- &
                  1._kp)/(alpha*(1._kp+1._kp/(2._kp*phi0**2))))

  end function hni_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function hni_efold_primitive(x,alpha,phi0)
    implicit none
    real(kp), intent(in) :: x,alpha,phi0
    real(kp) :: hni_efold_primitive

    hni_efold_primitive = -phi0**2/alpha*(log(tan(x/2._kp))+alpha*log(sin(x)))

  end function hni_efold_primitive



!returns x at bfold=-efolds before the end of inflation
  function hni_x_trajectory(bfold,xend,alpha,phi0)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,phi0
    real(kp) :: hni_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: hniData

    mini = epsilon(1._kp) !potential maximum
    maxi = xend*(1._kp-epsilon(1._kp))

    hniData%real1 = alpha
    hniData%real2 = phi0
    hniData%real3 = -bfold + hni_efold_primitive(xend,alpha,phi0)

    hni_x_trajectory = zbrent(find_hni_x_trajectory,mini,maxi,tolFind,hniData)

  end function hni_x_trajectory

  function find_hni_x_trajectory(x,hniData)
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: hniData
    real(kp) :: find_hni_x_trajectory
    real(kp) :: alpha,phi0,NplusPrimEnd

    alpha = hniData%real1
    phi0 = hniData%real2
    NplusPrimEnd = hniData%real3

    find_hni_x_trajectory = hni_efold_primitive(x,alpha,phi0) - NplusPrimEnd

  end function find_hni_x_trajectory



!Minimum value of alpha (depending on phi0) so that inflation ends naturally
  function hni_alphamin(phi0)
    implicit none
    real(kp) :: hni_alphamin
    real(kp), intent(in) :: phi0

    hni_alphamin = 1._kp/sqrt(1._kp+0.5_kp/phi0**2)

  end function hni_alphamin

!Maximum value of phi0 (depending on alpha) so that inflation ends naturally
  function hni_phizeromax(alpha)
    implicit none
    real(kp) :: hni_phizeromax
    real(kp), intent(in) :: alpha

    hni_phizeromax = alpha/(sqrt(2._kp*(1._kp - alpha**2)))

  end function hni_phizeromax



end module hnisr
