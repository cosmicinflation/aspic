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

module nmlfi1sr
  use infprec, only : kp,pi,tolkp,transfert
  use inftools, only : zbrent, easydverk
  use nmlfi1common, only : nmlfi1_x, nmlfi1_hbar, nmlfi1_deriv_x, nmlfi1_deriv_second_x
  use nmlfi1common, only : nmlfi1_norm_parametric_potential, nmlfi1_norm_deriv_parametric_potential
  use nmlfi1common, only : nmlfi1_norm_deriv_second_parametric_potential,nmlfi1_parametric_hbar_potmax
  use nmlfi1common, only : nmlfi1_parametric_epsilon_one, nmlfi1_parametric_epsilon_two
  use nmlfi1common, only : nmlfi1_parametric_epsilon_three, nmlfi1_parametric_efold_primitive
  use nmlfi1common, only : nmlfi1_parametric_hbar_trajectory, nmlfi1_parametric_hbar_endinf
  implicit none

  private

  public nmlfi1_norm_potential, nmlfi1_epsilon_one, nmlfi1_epsilon_two, nmlfi1_epsilon_three
  public nmlfi1_x_endinf, nmlfi1_efold_primitive, nmlfi1_x_trajectory
  public nmlfi1_norm_deriv_potential, nmlfi1_norm_deriv_second_potential 
  public nmlfi1_x, nmlfi1_hbar

contains
  


  function nmlfi1_x_endinf(xi,p)
    implicit none
    real(kp), intent(in) :: xi,p
    real(kp) :: nmlfi1_x_endinf

    real(kp) :: hbarend
    
    hbarend = nmlfi1_parametric_hbar_endinf(xi,p)

    nmlfi1_x_endinf = nmlfi_x(hbarend,xi)
  
  end function nmlfi1_x_endinf

 
!this is integral[V(phi)/V'(phi) dphi]
  function nmlfi1_efold_primitive(x,xi,p)
    implicit none
    real(kp), intent(in) :: x,xi,p
    real(kp) :: nmlfi1_efold_primitive
    
    nmlfi1_efold_primitive = nmlfi_efold_primitive(x,xi,p)

  end function nmlfi1_efold_primitive


  
  function nmlfi1_parametric_hbar_trajectory(bfold,hbarend,xi,p)
    implicit none
    real(kp), intent(in) :: bfold, hbarend, xi,p
    real(kp) :: nmlfi1_parametric_hbar_trajectory

    hbarmini = hbarend
    hbarmaxi = boo?
    
    nmlfi1_parametric_hbar_trajectory = nmlfi_parametric_hbar_trajectory(hbarend,hbarmini,hbarmaxi,xi,p)

  end function nmlfi1_parametric_hbar_trajectory

  
!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function nmlfi1_x_trajectory(bfold,xend,xi,p)
    implicit none
    real(kp), intent(in) :: bfold, xend, xi,p
    real(kp) :: nmlfi1_x_trajectory
    
    hbar = nmlfi1_parametric_hbar_trajectory(bfold,hbarend,xi,p)

    nmlfi1_x_trajectory = nmlfi_x(hbar,xi)
            
  end function nmlfi1_x_trajectory


end module nmlfi1sr
