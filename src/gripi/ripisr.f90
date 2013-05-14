!slow-roll functions for the renormalizable inflection point inflation potential
!
!V(phi) = M**4 (x**2 - 4/3 x**3 + 1/2 x**4)
!
!x = phi/phi0

module ripisr
  use infprec, only : kp,pi,tolkp,transfert
  use inftools, only : zbrent
  use gripicommon, only : gripi_norm_potential, gripi_norm_deriv_potential
  use gripicommon, only : gripi_norm_deriv_second_potential
  use gripicommon, only : gripi_epsilon_one, gripi_epsilon_two, gripi_epsilon_three
  use gripicommon, only : gripi_x_epsonemin, gripi_x_endinf

  implicit none
  
  private

  real(kp), parameter :: ripiAlpha = 1._kp

  public ripi_norm_potential, ripi_norm_deriv_potential
  public ripi_norm_deriv_second_potential
  public ripi_epsilon_one, ripi_epsilon_two, ripi_epsilon_three  
  public ripi_x_epsonemin, ripi_x_endinf

  public ripi_efold_primitive, ripi_x_trajectory
  

contains


  function ripi_norm_potential(x,phi0)
    implicit none
    real(kp) :: ripi_norm_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: phi0

    ripi_norm_potential = gripi_norm_potential(x,ripiAlpha,phi0)

  end function ripi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function ripi_norm_deriv_potential(x,phi0)
    implicit none
    real(kp) :: ripi_norm_deriv_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: phi0

    ripi_norm_deriv_potential = gripi_norm_deriv_potential(x,ripiAlpha,phi0)

  end function ripi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function ripi_norm_deriv_second_potential(x,phi0)
    implicit none
    real(kp) :: ripi_norm_deriv_second_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: phi0

    ripi_norm_deriv_second_potential = gripi_norm_deriv_second_potential(x,ripiAlpha,phi0)

  end function ripi_norm_deriv_second_potential



  !epsilon_one(x)
  function ripi_epsilon_one(x,phi0)    
    implicit none
    real(kp) :: ripi_epsilon_one
    real(kp), intent(in) :: x,phi0

    ripi_epsilon_one = gripi_epsilon_one(x,ripiAlpha,phi0)

  end function ripi_epsilon_one


  !epsilon_two(x)
  function ripi_epsilon_two(x,phi0)    
    implicit none
    real(kp) :: ripi_epsilon_two
    real(kp), intent(in) :: x,phi0

    ripi_epsilon_two = gripi_epsilon_two(x,ripiAlpha,phi0)

  end function ripi_epsilon_two


  !epsilon_three(x)
  function ripi_epsilon_three(x,phi0)    
    implicit none
    real(kp) :: ripi_epsilon_three
    real(kp), intent(in) :: x,phi0

    ripi_epsilon_three = gripi_epsilon_three(x,ripiAlpha,phi0)

  end function ripi_epsilon_three



  !returns x at the end of inflation defined as epsilon1=1
  function ripi_x_endinf(phi0)
    implicit none
    real(kp), intent(in) :: phi0
    real(kp) :: ripi_x_endinf

    ripi_x_endinf = gripi_x_endinf(ripiAlpha,phi0)

  end function ripi_x_endinf



  !this is integral(V(phi)/V'(phi) dphi)
  function ripi_efold_primitive(x,phi0)
    implicit none
    real(kp), intent(in) :: x,phi0
    real(kp) :: ripi_efold_primitive

    if (x .gt. 1._kp) then
       stop 'ripi_efold_primitive evaluated at x>=1 flat inflection point position!'
    endif

    ripi_efold_primitive = phi0**2*(-(2._kp*x/(9._kp)-x**2/8._kp+16._kp/(27._kp* &
         (3._kp*x-4._kp))+4._kp*log(4._kp-3._kp*x)/(27._kp)))

    ripi_efold_primitive = phi0**2 * ( &
         0.125_kp - 1._kp/(12._kp*(x-1._kp)) - x/6._kp + x**2/8._kp &
     - log(1._kp-x)/12._kp )
    

  end function ripi_efold_primitive


  !returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ripi_x_trajectory(bfold,xend,phi0)
    implicit none
    real(kp), intent(in) :: bfold, phi0, xend
    real(kp) :: ripi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ripiData


    mini = xEnd*(1._kp+epsilon(1._kp))
!Position of the flat inflection point
    maxi = 1._kp-epsilon(1._kp)

    ripiData%real1 = phi0
    ripiData%real2 = -bfold + ripi_efold_primitive(xend,phi0)

    ripi_x_trajectory = zbrent(find_ripi_x_trajectory,mini,maxi,tolFind,ripiData)

  end function ripi_x_trajectory

  function find_ripi_x_trajectory(x,ripiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ripiData
    real(kp) :: find_ripi_x_trajectory
    real(kp) :: phi0,NplusNuend

    phi0 = ripiData%real1
    NplusNuend = ripiData%real2

    find_ripi_x_trajectory = ripi_efold_primitive(x,phi0) - NplusNuend

  end function find_ripi_x_trajectory


  function ripi_x_epsonemin(phi0)
    implicit none
    real(kp) :: ripi_x_epsonemin, phi0

    ripi_x_epsonemin = gripi_x_epsonemin(ripiAlpha,phi0)

  end function ripi_x_epsonemin


end module ripisr
