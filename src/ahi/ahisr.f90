!slow-roll functions for the axion hilltop inflation potential
!
!V(phi) = M^4 * [ nu0 - 2 cos(phi/phi0) + (pi - phi/phi0) sin(phi/phi0) ]
!
!x = phi/phi0
!
!nu0 is a constant such that the potential vanishes at its two minimas around pi

module ahisr
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent, easydverk
  implicit none

  private

! Position of the lower mininum of the potential
  real(kp), parameter :: xmin = -1.3518168043192709368452375440008191980184144729_kp
! Position of the upper mininum of the potential
  real(kp), parameter :: xmax = 7.63500211149885741377052431055982496641275327166_kp
! Potential shift constant
  real(kp), parameter :: nu0 = 4.820572476962922009964861354665487247290442049596_kp

  public  ahi_norm_potential, ahi_epsilon_one, ahi_epsilon_two, ahi_epsilon_three
  public  ahi_x_endinf, ahi_efold_primitive, ahi_x_trajectory
  public  ahi_norm_deriv_potential, ahi_norm_deriv_second_potential


contains
!returns V/M^4
  function ahi_norm_potential(x, phi0)
    implicit none
    real(kp) :: ahi_norm_potential
    real(kp), intent(in) :: x, phi0

    ahi_norm_potential = nu0-2._kp*cos(x)+(acos(-1._kp)-x)*sin(x)

  end function ahi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function ahi_norm_deriv_potential(x, phi0)
    implicit none
    real(kp) :: ahi_norm_deriv_potential
    real(kp), intent(in) :: x, phi0

   ahi_norm_deriv_potential = ((acos(-1._kp)-x)*cos(x)+sin(x))

  end function ahi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function ahi_norm_deriv_second_potential(x, phi0)
    implicit none
    real(kp) :: ahi_norm_deriv_second_potential
    real(kp), intent(in) :: x, phi0

    ahi_norm_deriv_second_potential = (-acos(-1._kp)+x)*sin(x)

  end function ahi_norm_deriv_second_potential



!epsilon_one(x)
  function ahi_epsilon_one(x, phi0)    
    implicit none
    real(kp) :: ahi_epsilon_one
    real(kp), intent(in) :: x, phi0
    
    ahi_epsilon_one = ((acos(-1._kp)-x)*cos(x)+sin(x))**2/(2._kp*phi0**2* &
                    (nu0-2._kp*cos(x)+(acos(-1._kp)-x)*sin(x))**2)
    
  end function ahi_epsilon_one


!epsilon_two(x)
  function ahi_epsilon_two(x, phi0)    
    implicit none
    real(kp) :: ahi_epsilon_two
    real(kp), intent(in) :: x, phi0
    
    ahi_epsilon_two = (1._kp+2._kp*(acos(-1._kp)-x)**2-cos(2._kp*x)+ &
                      2._kp*nu0*(acos(-1._kp)-x)*sin(x))/(phi0**2* &
                      (nu0-2._kp*cos(x)+(acos(-1._kp)-x)*sin(x))**2)
    
  end function ahi_epsilon_two


!epsilon_three(x)
  function ahi_epsilon_three(x, phi0)    
    implicit none
    real(kp) :: ahi_epsilon_three
    real(kp), intent(in) :: x, phi0
    
    ahi_epsilon_three = (((acos(-1._kp)-x)*cos(x)+sin(x))*(9._kp*nu0* &
                        (acos(-1._kp)-x)+2*(-nu0**2+2*(-2._kp+(acos(-1._kp)-x)**2))* &
                        (acos(-1._kp)-x)*cos(x)+nu0*(-acos(-1._kp)+x)*cos(2*x)+ &
                        (5._kp+2._kp*nu0**2+8._kp*(acos(-1._kp)-x)**2)*sin(x)+ &
                        nu0*(-4._kp+(acos(-1._kp)-x)**2)*sin(2._kp*x)+sin(3._kp*x)))/ &
                        (phi0**2*(nu0-2._kp*cos(x)+(acos(-1._kp)-x)*sin(x))**2* &
                        (1._kp+2._kp*(acos(-1._kp)-x)**2-cos(2._kp*x)+2._kp*nu0* &
                        (acos(-1._kp)-x)*sin(x)))
    
  end function ahi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function ahi_x_endinf(phi0)
    implicit none
    real(kp), intent(in) :: phi0
    real(kp) :: ahi_x_endinf
    real(kp) :: mini, maxi
    real(kp), parameter :: tolFind=tolkp
    type(transfert) :: ahiData

    mini = xmin*(1._kp+epsilon(1._kp))
    maxi = acos(-1._kp)*(1._kp-epsilon(1._kp))

    ahiData%real1 = phi0

    ahi_x_endinf = zbrent(find_ahi_x_endinf,mini,maxi,tolFind,ahiData)
   
  end function ahi_x_endinf

  function find_ahi_x_endinf(x,ahiData)
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ahiData
    real(kp) :: find_ahi_x_endinf
    real(kp) :: phi0

    phi0 = ahiData%real1

    find_ahi_x_endinf = ahi_epsilon_one(x, phi0)-1._kp

  end function find_ahi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function ahi_efold_primitive(x, phi0)
    implicit none
    real(kp), intent(in) :: x, phi0
    real(kp) :: ahi_efold_primitive
    type(transfert) :: ahiData
    real(kp), parameter :: tolInt = tolkp
    integer, parameter :: neq = 1
    real(kp) :: xvar, xinf
    real(kp), dimension(neq) :: yvar

    !let us start where inflation ends
    xvar = ahi_x_endinf(phi0)
    yvar(1) = 0._kp

    ahiData%real1 = phi0

    call easydverk(neq,find_ahi_efold_primitive,xvar,yvar,x,tolInt,ahiData)

    ahi_efold_primitive = yvar(1)

  end function ahi_efold_primitive

  subroutine find_ahi_efold_primitive(n,x,y,yprime,ahiData)
    implicit none          
    integer :: n
    real(kp) :: x
    real(kp), dimension(n) :: y, yprime
    type(transfert), optional, intent(inout) :: ahiData
    real(kp) :: phi0

    phi0 = ahiData%real1

    yprime(1) = ahi_norm_potential(x, phi0)/ahi_norm_deriv_potential(x, phi0)*phi0**2

  end subroutine find_ahi_efold_primitive



!returns x at bfold=-efolds before the end of inflation
  function ahi_x_trajectory(bfold,xend,phi0)
    implicit none
    real(kp), intent(in) :: bfold,xend,phi0
    real(kp) :: ahi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ahiData

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = acos(-1._kp)*(1._kp-epsilon(1._kp))

    ahiData%real1 = phi0
    ahiData%real2 = -bfold + ahi_efold_primitive(xend,phi0)

    ahi_x_trajectory = zbrent(find_ahi_x_trajectory,mini,maxi,tolFind,ahiData)

  end function ahi_x_trajectory

  function find_ahi_x_trajectory(x,ahiData)
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ahiData
    real(kp) :: find_ahi_x_trajectory
    real(kp) :: phi0,NplusPrimEnd

    phi0 = ahiData%real1
    NplusPrimEnd = ahiData%real2

    find_ahi_x_trajectory = ahi_efold_primitive(x,phi0) - NplusPrimEnd

  end function find_ahi_x_trajectory



  
end module ahisr
