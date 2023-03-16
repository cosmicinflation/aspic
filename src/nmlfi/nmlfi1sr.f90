! Slow-roll functions for Non-Minimal Large Field Inflation 1 where
! inflation occurs at decreasing field values towards zero.
!
! This regime exists only provided inflation is supported, namely for
! p < p+ For 4 < p < p+, this is the only NMLFI regime possible. For p
! < 4, this regime occurs in the domain for which x < xVmax only.
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
  use infprec, only : kp
  use nmlficommon, only : hbarSmall, hbarBig, pplus
  use nmlficommon, only : nmlfi_norm_potential, nmlfi_norm_deriv_potential, nmlfi_norm_deriv_second_potential
  use nmlficommon, only : nmlfi_epsilon_one, nmlfi_epsilon_two, nmlfi_epsilon_three
  use nmlficommon, only : nmlfi_efold_primitive, nmlfi_parametric_hbar_trajectory
  use nmlficommon, only : nmlfi_hbarsquare_epsoneunity, nmlfi_hbar_potmax
  
  implicit none

  private

  public nmlfi1_norm_potential, nmlfi1_norm_deriv_potential, nmlfi1_norm_deriv_second_potential
  public nmlfi1_epsilon_one, nmlfi1_epsilon_two, nmlfi1_epsilon_three
  public nmlfi1_x_endinf, nmlfi1_efold_primitive, nmlfi1_x_trajectory
  public nmlfi1_hbar_endinf, nmlfi1_parametric_hbar_trajectory, nmlfi1_check_params


  
contains


  function nmlfi1_check_params(xi,p)
    implicit none
    logical :: nmlfi1_check_params
    real(kp), intent(in) :: xi,p

    nmlfi1_check_params = (p.lt.pplus)

  end function nmlfi1_check_params


  
  function nmlfi1_norm_potential(x,xi,p)
    implicit none
    real(kp) :: nmlfi1_norm_potential
    real(kp), intent(in) :: x,xi,p
    
    nmlfi1_norm_potential = nmlfi_norm_potential(hbar,xi,p)

  end function nmlfi1_norm_potential


  
!with respect to x
  function nmlfi1_norm_deriv_potential(x,xi,p)
    implicit none
    real(kp) :: nmlfi1_norm_deriv_potential
    real(kp), intent(in) :: x,xi,p

    nmlfi1_norm_deriv_potential = nmlfi_norm_deriv_potential(x,xi,p)

  end function nmlfi1_norm_deriv_potential


  !with respect to x
  function nmlfi1_norm_deriv_second_potential(x,xi,p)
    implicit none
    real(kp) :: nmlfi1_norm_deriv_second_potential
    real(kp), intent(in) :: x,xi,p
   
    nmlfi1_norm_deriv_second_potential = nmlfi_norm_deriv_second_potential(x,xi,p)

  end function nmlfi1_norm_deriv_second_potential

  
  function nmlfi_epsilon_one(x,xi,p)
    implicit none
    real(kp) :: nmlfi1_epsilon_one
    real(kp), intent(in) :: x,xi,p

    nmlfi1_epsilon_one = nmlfi_epsilon_one(x,xi,p)

  end function nmlfi_epsilon_one


  function nmlfi1_epsilon_two(x,xi,p)
    implicit none
    real(kp) :: nmlfi1_epsilon_two
    real(kp), intent(in) :: x,xi,p

    nmlfi1_epsilon_two = nmlfi_epsilon_two(x,xi,p)

  end function nmlfi1_epsilon_two


  function nmlfi1_epsilon_three(x,xi,p)
    implicit none
    real(kp) :: nmlfi1_epsilon_three
    real(kp), intent(in) :: x,xi,p

    nmlfi1_epsilon_three = nmlfi_epsilon_three(x,xi,p)

  end function nmlfi1_epsilon_three



  function nmlfi1_x_endinf(xi,p)
    implicit none
    real(kp) ::  nmlfi1_x_endinf

    real(kp) :: hbarendinf
    
    hbarendinf = nmlfi1_hbar_endinf(xi,p)

    nmlfi1_x_endinf = nmlfi_x(hbarendinf,xi)
    
  end function nmlfi1_x_endinf


  
  function nmlfi1_hbar_endinf(xi,p)
    implicit none
    real(kp) :: nmlfi1_hbar_endinf
    real(kp), intent(in) :: xi,p

    real(kp) :: hbarepsone2

    if (.not.nmlfi1_check_params(xi,p)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi1_hbar_endinf: NMLFI1 does not exist for these parameters!'
    endif
    
    hbarepsone2 = nmlfi_hbarsquare_epsoneunity(xi,p)

!in all cases, the hbarepsone2(1) is the end of inflation at x < xVmax
!if xVmax exists

!safety check    
    if (hbarepsone2(1).le.0._kp) stop 'nmlfi1_hbar_endinf: hbarend^2 < 0!?'

    nmlfi1_hbar_endinf = sqrt(hbarepsone2(1))
        
  end function nmlfi1_hbar_endinf

  
  
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
    real(kp) :: hbarmin, hbarmax

    if (.not.nmlfi1_check_params(xi,p)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi1_parametric_hbar_trajectory: NMLFI1 does not exist for these parameters!'
    endif
    
    hbarmin = hbarend

    if (p.ge.4._kp) then
       hbarmax = hbarBig
    else
       hbarmax = nmlfi_hbar_potmax(xi,p)
    end if
    
    nmlfi1_parametric_hbar_trajectory = nmlfi_parametric_hbar_trajectory(bfold,hbarend,hbarmin,hbarmax,xi,p)

  end function nmlfi1_parametric_hbar_trajectory

  
!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function nmlfi1_x_trajectory(bfold,xend,xi,p)
    implicit none
    real(kp), intent(in) :: bfold, xend, xi,p
    real(kp) :: nmlfi1_x_trajectory
    real(kp) :: hbar
    
    hbar = nmlfi1_parametric_hbar_trajectory(bfold,hbarend,xi,p)

    nmlfi1_x_trajectory = nmlfi_x(hbar,xi)
            
  end function nmlfi1_x_trajectory


end module nmlfi1sr
