! Slow-roll functions for Non-Minimal Large Field Inflation
!
! V(hbar) = M^4[ hbar^p / (1+hbar^2)^2 ]
!
! x = phi/Mg ~ phi/Mpl, the Einstein frame field value in unit of the
!
! Jordan Frame gravity scale Mg
!
! (dx/dhbar)^2 = [ 1 + (1+6xi) hbar^2 ] / [ xi (1+hbar^2)^2]
!
! NB: hbar = sqrt(xi) phibar/Mg where phibar is the Jordan frame field
!
! The normalization of the potential is:
!
!  M^4 = lambdabar Mg^4 / xi^2
!
! all with Mg ~ Mpl.

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
  function nmlfi_norm_potential(x,xi,p)
    implicit none
    real(kp) :: nmlfi_norm_potential
    real(kp), intent(in) :: x,xi,p
    real(kp) :: hbar    

    hbar = nmlfi_hbar(x,xi)
    
    nmlfi_norm_potential = nmlfi_norm_parametric_potential(hbar,xi,p)

  end function nmlfi_norm_potential



!with respect to x
  function nmlfi_norm_deriv_potential(x,xi,p)
    implicit none
    real(kp) :: nmlfi_norm_deriv_potential
    real(kp), intent(in) :: x,xi,p
    real(kp) :: hbar, dx

    hbar = nmlfi_hbar(x,xi)
    
    dx = nmlfi_deriv_x(hbar,xi)

    nmlfi_norm_deriv_potential = nmlfi_norm_deriv_parametric_potential(hbar,xi,p)/dx

  end function nmlfi_norm_deriv_potential



!with respect to x
  function nmlfi_norm_deriv_second_potential(x,xi,p)
    implicit none
    real(kp) :: nmlfi_norm_deriv_second_potential
    real(kp), intent(in) :: x,xi,p
    real(kp) :: hbar, dV, d2V, dx, d2x

    hbar = nmlfi_hbar(x,xi)
    
    dV = nmlfi_norm_deriv_parametric_potential(hbar,xi,p)
    d2V = nmlfi_norm_deriv_second_parametric_potential(hbar,xi,p)

    dx = nmlfi_deriv_x(hbar,xi)
    d2x = nmlfi_deriv_second_x(hbar,xi)

    nmlfi_norm_deriv_second_potential = (d2V - (d2x/dx)*dV)/dx/dx

  end function nmlfi_norm_deriv_second_potential


  
  function nmlfi_epsilon_one(x,xi,p)
    implicit none
    real(kp) :: nmlfi_epsilon_one
    real(kp), intent(in) :: x,xi,p
    real(kp) :: hbar

    hbar = nmlfi_hbar(x,xi)

    nmlfi_epsilon_one = nmlfi_parametric_epsilon_one(hbar,xi,p)

  end function nmlfi_epsilon_one



  function nmlfi_epsilon_two(x,xi,p)
    implicit none
    real(kp) :: nmlfi_epsilon_two
    real(kp), intent(in) :: x,xi,p
    real(kp) :: hbar

    hbar = nmlfi_hbar(x,xi)

    nmlfi_epsilon_two = nmlfi_parametric_epsilon_two(hbar,xi,p)

  end function nmlfi_epsilon_two


  function nmlfi_epsilon_three(x,xi,p)
    implicit none
    real(kp) :: nmlfi_epsilon_three
    real(kp), intent(in) :: x,xi,p

    real(kp) :: hbar

    hbar = nmlfi_hbar(x,xi)

    nmlfi_epsilon_three = nmlfi_parametric_epsilon_three(hbar,xi,p)

  end function nmlfi_epsilon_three



  function nmlfi_x_endinf(xi,p)
    implicit none
    real(kp), intent(in) :: xi,p
    real(kp) :: nmlfi_x_endinf

    real(kp) :: hbarend
    
    hbarend = nmlfi_parametric_hbar_endinf(xi,p)

    nmlfi_x_endinf = nmlfi_x(hbarend,xi)
  
  end function nmlfi_x_endinf


!the field value at which the potential is maximal (exists only for p<4)
  function nmlfi_x_potmax(xi,p)
    implicit none
    real(kp), intent(in) :: xi,p
    real(kp) :: nmlfi_x_potmax

    real(kp) :: hbarmax

    if (p.ge.4._kp) then
       stop 'nmlfi_x_potmax: p>=4 no maximum'
    endif
    
    hbarmax = nmlfi_parametric_hbar_potmax(xi,p)

    nmlfi_x_potmax = nmlfi_x(hbarmax,xi)
  
  end function nmlfi_x_potmax
  

 
!this is integral[V(phi)/V'(phi) dphi]
  function nmlfi_efold_primitive(x,xi,p)
    implicit none
    real(kp), intent(in) :: x,xi,p
    real(kp) :: nmlfi_efold_primitive

    real(kp) :: hbar

    hbar = nmlfi_hbar(x,xi)

    nmlfi_efold_primitive = nmlfi_parametric_efold_primitive(hbar,xi,p)

  end function nmlfi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function nmlfi_x_trajectory(bfold,xend,xi,p)
    implicit none
    real(kp), intent(in) :: bfold, xend, xi,p
    real(kp) :: nmlfi_x_trajectory
    real(kp) :: hbar,hbarend

    hbarend = nmlfi_parametric_hbar_endinf(xi,p)

    hbar = nmlfi_parametric_hbar_trajectory(bfold,hbarend,xi,p)

    nmlfi_x_trajectory = nmlfi_x(hbar,xi)

  end function nmlfi_x_trajectory


end module nmlfisr
