! Slow-roll functions for Higgs inflation out of the large field approximation
!
! V(hbar) = M^4[ (hbar^2 - xi v^2)^2 / (1+hbar^2)^2 ]
!
! x = phi/Mg ~ phi/Mpl, the field value in unit of the Jordan Frame gravity scale Mg
!
! (dx/dhbar)^2 = [ 1 + (1+6xi) hbar^2 ] / [ xi (1+hbar^2)^2]
!
! NB: hbar = sqrt(xi) H/Mg where H is the Higgs field in the unitary gauge
!
! The normalization of the potential is:
!
!  M^4 = lambda Mg^4 / (4 xi^2)
!
! with Mg ~ Mpl.

module hisr
  use infprec, only : kp,pi,tolkp,transfert
  use inftools, only : zbrent, easydverk
  use hicommon, only : hi_x, hi_hbar, hi_deriv_x, hi_deriv_second_x
  use hicommon, only : hi_norm_parametric_potential, hi_norm_deriv_parametric_potential
  use hicommon, only : hi_norm_deriv_second_parametric_potential
  use hicommon, only : hi_parametric_epsilon_one, hi_parametric_epsilon_two
  use hicommon, only : hi_parametric_epsilon_three, hi_parametric_efold_primitive
  use hicommon, only : hi_parametric_hbar_trajectory, hi_parametric_hbar_endinf
  implicit none

  private

  public hi_norm_potential, hi_epsilon_one, hi_epsilon_two, hi_epsilon_three
  public hi_x_endinf, hi_efold_primitive, hi_x_trajectory
  public hi_norm_deriv_potential, hi_norm_deriv_second_potential 
  public hi_x, hi_hbar

contains
  


!returns V/M^4
  function hi_norm_potential(x,xi)
    implicit none
    real(kp) :: hi_norm_potential
    real(kp), intent(in) :: x,xi
    real(kp) :: hbar    

    hbar = hi_hbar(x,xi)
    
    hi_norm_potential = hi_norm_parametric_potential(hbar,xi)

  end function hi_norm_potential



!with respect to x
  function hi_norm_deriv_potential(x,xi)
    implicit none
    real(kp) :: hi_norm_deriv_potential
    real(kp), intent(in) :: x,xi
    real(kp) :: hbar, dx

    hbar = hi_hbar(x,xi)
    
    dx = hi_deriv_x(hbar,xi)

    hi_norm_deriv_potential = hi_norm_deriv_parametric_potential(hbar,xi)/dx

  end function hi_norm_deriv_potential



!with respect to x
  function hi_norm_deriv_second_potential(x,xi)
    implicit none
    real(kp) :: hi_norm_deriv_second_potential
    real(kp), intent(in) :: x,xi
    real(kp) :: hbar, dV, d2V, dx, d2x

    hbar = hi_hbar(x,xi)
    
    dV = hi_norm_deriv_parametric_potential(hbar,xi)
    d2V = hi_norm_deriv_second_parametric_potential(hbar,xi)

    dx = hi_deriv_x(hbar,xi)
    d2x = hi_deriv_second_x(hbar,xi)

    hi_norm_deriv_second_potential = (d2V - (d2x/dx)*dV)/dx/dx

  end function hi_norm_deriv_second_potential


  
  function hi_epsilon_one(x,xi)
    implicit none
    real(kp) :: hi_epsilon_one
    real(kp), intent(in) :: x,xi
    real(kp) :: hbar

    hbar = hi_hbar(x,xi)

    hi_epsilon_one = hi_parametric_epsilon_one(hbar,xi)

  end function hi_epsilon_one



  function hi_epsilon_two(x,xi)
    implicit none
    real(kp) :: hi_epsilon_two
    real(kp), intent(in) :: x,xi
    real(kp) :: hbar

    hbar = hi_hbar(x,xi)

    hi_epsilon_two = hi_parametric_epsilon_two(hbar,xi)

  end function hi_epsilon_two


  function hi_epsilon_three(x,xi)
    implicit none
    real(kp) :: hi_epsilon_three
    real(kp), intent(in) :: x,xi

    real(kp) :: hbar

    hbar = hi_hbar(x,xi)

    hi_epsilon_three = hi_parametric_epsilon_three(hbar,xi)

  end function hi_epsilon_three



  function hi_x_endinf(xi)
    implicit none
    real(kp), intent(in) :: xi
    real(kp) :: hi_x_endinf

    real(kp) :: hbarend
    
    hbarend = hi_parametric_hbar_endinf(xi)

    hi_x_endinf = hi_x(hbarend,xi)
  
  end function hi_x_endinf


 
!this is integral[V(phi)/V'(phi) dphi]
  function hi_efold_primitive(x,xi)
    implicit none
    real(kp), intent(in) :: x,xi
    real(kp) :: hi_efold_primitive

    real(kp) :: hbar

    hbar = hi_hbar(x,xi)

    hi_efold_primitive = hi_parametric_efold_primitive(hbar,xi)

  end function hi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function hi_x_trajectory(bfold,xend,xi)
    implicit none
    real(kp), intent(in) :: bfold, xend, xi
    real(kp) :: hi_x_trajectory
    real(kp) :: hbar,hbarend

    hbarend = hi_parametric_hbar_endinf(xi)

    hbar = hi_parametric_hbar_trajectory(bfold,hbarend,xi)

    hi_x_trajectory = hi_x(hbar,xi)

  end function hi_x_trajectory


end module hisr
