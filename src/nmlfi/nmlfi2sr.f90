! Slow-roll functions for Non-Minimal Large Field Inflation 2 where
! inflation occurs at increasing field values for x>xVmax and
! gracefully ends at epsilon=1.
!
! This regime exists only under the conditions: [ p < p_ and xi > xizero(p) ]
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

module nmlfi2sr
  use infprec, only : kp
  use nmlficommon, only : pminus, pplus, nmlfi_hbar, nmlfi_x, nmlfi_xizero
  use nmlficommon, only : nmlfi_norm_potential, nmlfi_norm_deriv_potential, nmlfi_norm_deriv_second_potential
  use nmlficommon, only : nmlfi_epsilon_one, nmlfi_epsilon_two, nmlfi_epsilon_three
  use nmlficommon, only : nmlfi_efold_primitive, nmlfi_parametric_hbar_trajectory
  use nmlficommon, only : nmlfi_hbarsquare_epsoneunity, nmlfi_hbar_potmax
  use nmlficommon, only : nmlfi_numacc_hbarsquare_epsonenull, nmlfi_parametric_efold_primitive
  
  implicit none

  private

  public nmlfi2_check_params
  public nmlfi2_norm_potential, nmlfi2_norm_deriv_potential, nmlfi2_norm_deriv_second_potential
  public nmlfi2_epsilon_one, nmlfi2_epsilon_two, nmlfi2_epsilon_three
  public nmlfi2_x_endinf, nmlfi2_efold_primitive, nmlfi2_x_trajectory
  public nmlfi2_hbar_endinf, nmlfi2_parametric_hbar_trajectory

  public nmlfi2_numacc_hbarinimin, nmlfi2_numacc_xinimin, nmlfi2_numacc_efoldmax

  
contains


  function nmlfi2_check_params(xi,p)
    implicit none
    logical :: nmlfi2_check_params
    real(kp), intent(in) :: xi,p
    
    nmlfi2_check_params = (p.lt.pminus).and.(xi.ge.nmlfi_xizero(p))

  end function nmlfi2_check_params


  
  function nmlfi2_norm_potential(x,xi,p)
    implicit none
    real(kp) :: nmlfi2_norm_potential
    real(kp), intent(in) :: x,xi,p
    
    nmlfi2_norm_potential = nmlfi_norm_potential(x,xi,p)

  end function nmlfi2_norm_potential


  
!with respect to x
  function nmlfi2_norm_deriv_potential(x,xi,p)
    implicit none
    real(kp) :: nmlfi2_norm_deriv_potential
    real(kp), intent(in) :: x,xi,p

    nmlfi2_norm_deriv_potential = nmlfi_norm_deriv_potential(x,xi,p)

  end function nmlfi2_norm_deriv_potential


  !with respect to x
  function nmlfi2_norm_deriv_second_potential(x,xi,p)
    implicit none
    real(kp) :: nmlfi2_norm_deriv_second_potential
    real(kp), intent(in) :: x,xi,p
   
    nmlfi2_norm_deriv_second_potential = nmlfi_norm_deriv_second_potential(x,xi,p)

  end function nmlfi2_norm_deriv_second_potential

  
  function nmlfi2_epsilon_one(x,xi,p)
    implicit none
    real(kp) :: nmlfi2_epsilon_one
    real(kp), intent(in) :: x,xi,p

    nmlfi2_epsilon_one = nmlfi_epsilon_one(x,xi,p)

  end function nmlfi2_epsilon_one


  function nmlfi2_epsilon_two(x,xi,p)
    implicit none
    real(kp) :: nmlfi2_epsilon_two
    real(kp), intent(in) :: x,xi,p

    nmlfi2_epsilon_two = nmlfi_epsilon_two(x,xi,p)

  end function nmlfi2_epsilon_two


  function nmlfi2_epsilon_three(x,xi,p)
    implicit none
    real(kp) :: nmlfi2_epsilon_three
    real(kp), intent(in) :: x,xi,p

    nmlfi2_epsilon_three = nmlfi_epsilon_three(x,xi,p)

  end function nmlfi2_epsilon_three



  function nmlfi2_x_endinf(xi,p)
    implicit none
    real(kp) ::  nmlfi2_x_endinf
    real(kp), intent(in) :: xi,p
    
    real(kp) :: hbarendinf

    hbarendinf = nmlfi2_hbar_endinf(xi,p)

    nmlfi2_x_endinf = nmlfi_x(hbarendinf,xi)
    
  end function nmlfi2_x_endinf


  
  function nmlfi2_hbar_endinf(xi,p)
    implicit none
    real(kp) :: nmlfi2_hbar_endinf
    real(kp), intent(in) :: xi,p

    real(kp), dimension(2) :: hbarepsone2

    if (.not.nmlfi2_check_params(xi,p)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi2_hbar_endinf: NMLFI2 does not exist for these parameters!'
    endif
    
    hbarepsone2 = nmlfi_hbarsquare_epsoneunity(xi,p)

!safety check    
    if (hbarepsone2(2).le.0._kp) then
       write(*,*)'xi= p= hbarepsone2= ',xi,p,hbarepsone2
       stop 'nmlfi2_hbar_endinf: hbarend^2 < 0!?'
    endif
    
    nmlfi2_hbar_endinf = sqrt(hbarepsone2(2))
        
  end function nmlfi2_hbar_endinf

  


!A the top of the potential, eps1->0. This function returns the
!minimal possible value of hbarini, close to the top of the potential,
!starting at eps1=machine precision
  function nmlfi2_numacc_hbarinimin(xi,p)
    implicit none
    real(kp) :: nmlfi2_numacc_hbarinimin
    real(kp), intent(in) :: xi,p

    real(kp), dimension(2) :: hbar2null

    if (.not.nmlfi2_check_params(xi,p)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi2_numacc_hbarinimin: NMLFI2 does not exist for these parameters!'
    endif
    
    hbar2null = nmlfi_numacc_hbarsquare_epsonenull(xi,p)

    if (hbar2null(2).le.0._kp) then
       write(*,*)'xi= p= hbar2null= ',xi,p,hbar2null(:)
       stop ' nmlfi2_numacc_hbarinimin: bug!?'
    endif

    nmlfi2_numacc_hbarinimin = sqrt(hbar2null(2))

!safety check
    if (nmlfi2_numacc_hbarinimin.lt.nmlfi_hbar_potmax(xi,p)) then
       stop 'nmlfi2_numacc_hbarinimin: hbarinimin <  hbarpotmax!?'
    endif

  end function nmlfi2_numacc_hbarinimin
    
  

  function nmlfi2_numacc_xinimin(xi,p)
    implicit none
    real(kp) :: nmlfi2_numacc_xinimin
    real(kp), intent(in) :: xi,p
    real(kp) :: hbarinimin

    hbarinimin = nmlfi2_numacc_hbarinimin(xi,p)

    nmlfi2_numacc_xinimin = nmlfi_x(hbarinimin,xi)
    
  end function nmlfi2_numacc_xinimin


!the maximal number of efold computable at the current numerical accuracy
  function nmlfi2_numacc_efoldmax(xi,p)
    implicit none
    real(kp) :: nmlfi2_numacc_efoldmax
    real(kp), intent(in) :: xi,p

    real(kp) :: hbarend, hbarinimin
    
    hbarend = nmlfi2_hbar_endinf(xi,p)
    hbarinimin = nmlfi2_numacc_hbarinimin(xi,p)
    
    nmlfi2_numacc_efoldmax = -nmlfi_parametric_efold_primitive(hbarend,xi,p) &
         + nmlfi_parametric_efold_primitive(hbarinimin,xi,p)
    
  end function nmlfi2_numacc_efoldmax

  
  
  !this is integral[V(phi)/V'(phi) dphi]
  function nmlfi2_efold_primitive(x,xi,p)
    implicit none
    real(kp), intent(in) :: x,xi,p
    real(kp) :: nmlfi2_efold_primitive

    nmlfi2_efold_primitive = nmlfi_efold_primitive(x,xi,p)

  end function nmlfi2_efold_primitive



  function nmlfi2_parametric_hbar_trajectory(bfold,hbarend,xi,p)
    implicit none
    real(kp), intent(in) :: bfold, hbarend, xi,p
    real(kp) :: nmlfi2_parametric_hbar_trajectory
    real(kp) :: hbarmin, hbarmax

    if (.not.nmlfi2_check_params(xi,p)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi2_parametric_hbar_trajectory: NMLFI2 does not exist for these parameters!'
    endif
       
    hbarmin = nmlfi_hbar_potmax(xi,p)

    hbarmax = nmlfi2_hbar_endinf(xi,p)
        
    nmlfi2_parametric_hbar_trajectory = nmlfi_parametric_hbar_trajectory(bfold,hbarend,hbarmin,hbarmax,xi,p)

  end function nmlfi2_parametric_hbar_trajectory

  
!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function nmlfi2_x_trajectory(bfold,xend,xi,p)
    implicit none
    real(kp), intent(in) :: bfold, xend, xi,p
    real(kp) :: nmlfi2_x_trajectory
    real(kp) :: hbar, hbarend

    hbarend = nmlfi_hbar(xend,xi)
    
    hbar = nmlfi2_parametric_hbar_trajectory(bfold,hbarend,xi,p)

    nmlfi2_x_trajectory = nmlfi_x(hbar,xi)
            
  end function nmlfi2_x_trajectory


end module nmlfi2sr
