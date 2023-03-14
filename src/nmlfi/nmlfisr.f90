! Slow-roll functions for Higgs inflation out of the large field approximation
!
! V(hbar) = M^4[ (hbar^2 - xi v^2)^2 / (1+hbar^2) ]
!
! x = pnmlfi/Mg ~ pnmlfi/Mpl, the field value in unit of the Jordan Frame gravity scale Mg
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

module nmlfisr
  use infprec, only : kp,pi,tolkp,transfert
  use inftools, only : zbrent, easydverk
  use nmlficommon, only : nmlfi_x, nmlfi_hbar, nmlfi_deriv_x, nmlfi_deriv_second_x
  use nmlficommon, only : nmlfi_norm_parametric_potential, nmlfi_norm_deriv_parametric_potential
  use nmlficommon, only : nmlfi_norm_deriv_second_parametric_potential
  use nmlficommon, only : nmlfi_parametric_epsilon_one, nmlfi_parametric_epsilon_two
  use nmlficommon, only : nmlfi_parametric_epsilon_three, nmlfi_parametric_efold_primitive
  use nmlficommon, only : nmlfi_parametric_hbar_trajectory, nmlfi_parametric_hbar_endinf
  implicit none

  private

  public nmlfi_norm_potential, nmlfi_epsilon_one, nmlfi_epsilon_two, nmlfi_epsilon_three
  public nmlfi_x_endinf, nmlfi_efold_primitive, nmlfi_x_trajectory
  public nmlfi_norm_deriv_potential, nmlfi_norm_deriv_second_potential 
  public nmlfi_x, nmlfi_hbar

contains
  


!returns V/M^4
  function nmlfi_norm_potential(x,xi)
    implicit none
    real(kp) :: nmlfi_norm_potential
    real(kp), intent(in) :: x,xi
    real(kp) :: hbar    

    hbar = nmlfi_hbar(x,xi)
    
    nmlfi_norm_potential = nmlfi_norm_parametric_potential(hbar,xi)

  end function nmlfi_norm_potential



!with respect to x
  function nmlfi_norm_deriv_potential(x,xi)
    implicit none
    real(kp) :: nmlfi_norm_deriv_potential
    real(kp), intent(in) :: x,xi
    real(kp) :: hbar, dx

    hbar = nmlfi_hbar(x,xi)
    
    dx = nmlfi_deriv_x(hbar,xi)

    nmlfi_norm_deriv_potential = nmlfi_norm_deriv_parametric_potential(hbar,xi)/dx

  end function nmlfi_norm_deriv_potential



!with respect to x
  function nmlfi_norm_deriv_second_potential(x,xi)
    implicit none
    real(kp) :: nmlfi_norm_deriv_second_potential
    real(kp), intent(in) :: x,xi
    real(kp) :: hbar, dV, d2V, dx, d2x

    hbar = nmlfi_hbar(x,xi)
    
    dV = nmlfi_norm_deriv_parametric_potential(hbar,xi)
    d2V = nmlfi_norm_deriv_second_parametric_potential(hbar,xi)

    dx = nmlfi_deriv_x(hbar,xi)
    d2x = nmlfi_deriv_second_x(hbar,xi)

    nmlfi_norm_deriv_second_potential = (d2V - (d2x/dx)*dV)/dx/dx

  end function nmlfi_norm_deriv_second_potential


  
  function nmlfi_epsilon_one(x,xi)
    implicit none
    real(kp) :: nmlfi_epsilon_one
    real(kp), intent(in) :: x,xi
    real(kp) :: hbar

    hbar = nmlfi_hbar(x,xi)

    nmlfi_epsilon_one = nmlfi_parametric_epsilon_one(hbar,xi)

  end function nmlfi_epsilon_one



  function nmlfi_epsilon_two(x,xi)
    implicit none
    real(kp) :: nmlfi_epsilon_two
    real(kp), intent(in) :: x,xi
    real(kp) :: hbar

    hbar = nmlfi_hbar(x,xi)

    nmlfi_epsilon_two = nmlfi_parametric_epsilon_two(hbar,xi)

  end function nmlfi_epsilon_two


  function nmlfi_epsilon_three(x,xi)
    implicit none
    real(kp) :: nmlfi_epsilon_three
    real(kp), intent(in) :: x,xi

    real(kp) :: hbar

    hbar = nmlfi_hbar(x,xi)

    nmlfi_epsilon_three = nmlfi_parametric_epsilon_three(hbar,xi)

  end function nmlfi_epsilon_three



  function nmlfi_x_endinf(xi)
    implicit none
    real(kp), intent(in) :: xi
    real(kp) :: nmlfi_x_endinf

    real(kp) :: hbarend
    
    hbarend = nmlfi_parametric_hbar_endinf(xi)

    nmlfi_x_endinf = nmlfi_x(hbarend,xi)
  
  end function nmlfi_x_endinf


 
!tnmlfis is integral[V(pnmlfi)/V'(pnmlfi) dpnmlfi]
  function nmlfi_efold_primitive(x,xi)
    implicit none
    real(kp), intent(in) :: x,xi
    real(kp) :: nmlfi_efold_primitive

    real(kp) :: hbar

    hbar = nmlfi_hbar(x,xi)

    nmlfi_efold_primitive = nmlfi_parametric_efold_primitive(hbar,xi)

  end function nmlfi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function nmlfi_x_trajectory(bfold,xend,xi)
    implicit none
    real(kp), intent(in) :: bfold, xend, xi
    real(kp) :: nmlfi_x_trajectory
    real(kp) :: hbar,hbarend

    hbarend = nmlfi_parametric_hbar_endinf(xi)

    hbar = nmlfi_parametric_hbar_trajectory(bfold,hbarend,xi)

    nmlfi_x_trajectory = nmlfi_x(hbar,xi)

  end function nmlfi_x_trajectory


end module nmlfisr
