! Slow-roll functions for Non-Minimal Large Field Inflation 3 where
! inflation occurs at increasing field values for x>xVmax and does not
! gracefully ends. xend is an extra model parameter
!
! This regime exists only under two conditions: [ p < p_ and xi < xizero(p) ] or [ p_ < p < 4]
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

module nmlfi3sr
  use infprec, only : kp
  use nmlficommon, only : hbarBig, pminus, nmlfi_xizero, nmlfi_hbar, nmlfi_x
  use nmlficommon, only : nmlfi_norm_potential, nmlfi_norm_deriv_potential, nmlfi_norm_deriv_second_potential
  use nmlficommon, only : nmlfi_epsilon_one, nmlfi_epsilon_two, nmlfi_epsilon_three
  use nmlficommon, only : nmlfi_efold_primitive, nmlfi_parametric_hbar_trajectory
  use nmlficommon, only : nmlfi_hbar_potmax, nmlfi_numacc_hbarsquare_epsonenull
  use nmlficommon, only : nmlfi_parametric_efold_primitive
  
  implicit none

  private

  public nmlfi3_check_params
  public nmlfi3_norm_potential, nmlfi3_norm_deriv_potential, nmlfi3_norm_deriv_second_potential
  public nmlfi3_epsilon_one, nmlfi3_epsilon_two, nmlfi3_epsilon_three
  public nmlfi3_efold_primitive, nmlfi3_x_trajectory
  public nmlfi3_parametric_hbar_trajectory

  public nmlfi3_numacc_hbarinimin, nmlfi3_numacc_xinimin, nmlfi3_numacc_efoldmax
  public nmlfi3_numacc_hbarendmin
  
contains


  function nmlfi3_check_params(xi,p)
    implicit none
    logical :: nmlfi3_check_params
    real(kp), intent(in) :: xi,p

    nmlfi3_check_params = ( (p.le.pminus).and.(xi.lt.nmlfi_xizero(p)) ) &
         .or. ( (p.gt.pminus).and.(p.lt.4._kp) ) 

  end function nmlfi3_check_params


  
  function nmlfi3_norm_potential(x,xi,p)
    implicit none
    real(kp) :: nmlfi3_norm_potential
    real(kp), intent(in) :: x,xi,p
    
    nmlfi3_norm_potential = nmlfi_norm_potential(x,xi,p)

  end function nmlfi3_norm_potential


  
!with respect to x
  function nmlfi3_norm_deriv_potential(x,xi,p)
    implicit none
    real(kp) :: nmlfi3_norm_deriv_potential
    real(kp), intent(in) :: x,xi,p

    nmlfi3_norm_deriv_potential = nmlfi_norm_deriv_potential(x,xi,p)

  end function nmlfi3_norm_deriv_potential


  !with respect to x
  function nmlfi3_norm_deriv_second_potential(x,xi,p)
    implicit none
    real(kp) :: nmlfi3_norm_deriv_second_potential
    real(kp), intent(in) :: x,xi,p
   
    nmlfi3_norm_deriv_second_potential = nmlfi_norm_deriv_second_potential(x,xi,p)

  end function nmlfi3_norm_deriv_second_potential

  
  function nmlfi3_epsilon_one(x,xi,p)
    implicit none
    real(kp) :: nmlfi3_epsilon_one
    real(kp), intent(in) :: x,xi,p

    nmlfi3_epsilon_one = nmlfi_epsilon_one(x,xi,p)

  end function nmlfi3_epsilon_one


  function nmlfi3_epsilon_two(x,xi,p)
    implicit none
    real(kp) :: nmlfi3_epsilon_two
    real(kp), intent(in) :: x,xi,p

    nmlfi3_epsilon_two = nmlfi_epsilon_two(x,xi,p)

  end function nmlfi3_epsilon_two


  function nmlfi3_epsilon_three(x,xi,p)
    implicit none
    real(kp) :: nmlfi3_epsilon_three
    real(kp), intent(in) :: x,xi,p

    nmlfi3_epsilon_three = nmlfi_epsilon_three(x,xi,p)

  end function nmlfi3_epsilon_three


  
!A the top of the potential, eps1->0. This function returns the
!minimal possible value of hbarini, close to the top of the potential,
!starting at eps1=machine precision
  function nmlfi3_numacc_hbarinimin(xi,p)
    implicit none
    real(kp) :: nmlfi3_numacc_hbarinimin
    real(kp), intent(in) :: xi,p

    real(kp), dimension(2) :: hbar2null

    if (.not.nmlfi3_check_params(xi,p)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi3_numacc_hbarinimin: NMLFI3 does not exist for these parameters!'
    endif
    
    hbar2null = nmlfi_numacc_hbarsquare_epsonenull(xi,p)

    if (hbar2null(2).le.0._kp) then
       write(*,*)'xi= p= hbar2null= ',xi,p,hbar2null(:)
       stop ' nmlfi3_numacc_hbarinimin: bug!?'
    endif

    nmlfi3_numacc_hbarinimin = sqrt(hbar2null(2))

!safety check
    if (nmlfi3_numacc_hbarinimin.lt.nmlfi_hbar_potmax(xi,p)) then
       stop 'nmlfi3_numacc_hbarinimin: hbarinimin <  hbarpotmax!?'
    endif

  end function nmlfi3_numacc_hbarinimin
    
  

  function nmlfi3_numacc_xinimin(xi,p)
    implicit none
    real(kp) :: nmlfi3_numacc_xinimin
    real(kp), intent(in) :: xi,p
    real(kp) :: hbarinimin

    hbarinimin = nmlfi3_numacc_hbarinimin(xi,p)

    nmlfi3_numacc_xinimin = nmlfi_x(hbarinimin,xi)
    
  end function nmlfi3_numacc_xinimin


  function nmlfi3_numacc_hbarendmin(efold,xi,p)
    implicit none
    real(kp) :: nmlfi3_numacc_hbarendmin
    real(kp), intent(in) :: efold,xi,p

    real(kp) :: hbarinimin

    hbarinimin = nmlfi3_numacc_hbarinimin(xi,p)

    nmlfi3_numacc_hbarendmin = nmlfi3_parametric_hbar_trajectory(efold,hbarinimin,xi,p)

  end function nmlfi3_numacc_hbarendmin

    
  
!the maximal number of efold computable at the current numerical
!accuracy. This is not necessarily huge as eps1 can be not small
!asymptotically
  function nmlfi3_numacc_efoldmax(xi,p)
    implicit none
    real(kp) :: nmlfi3_numacc_efoldmax
    real(kp), intent(in) :: xi,p

    real(kp) :: hbarinimin, hbarendmax
    
    hbarendmax = huge(1._kp)
    hbarinimin = nmlfi3_numacc_hbarinimin(xi,p)
    
    nmlfi3_numacc_efoldmax = -nmlfi_parametric_efold_primitive(hbarendmax,xi,p) &
         + nmlfi_parametric_efold_primitive(hbarinimin,xi,p)
    
  end function nmlfi3_numacc_efoldmax


    
  
  !this is integral[V(phi)/V'(phi) dphi]
  function nmlfi3_efold_primitive(x,xi,p)
    implicit none
    real(kp), intent(in) :: x,xi,p
    real(kp) :: nmlfi3_efold_primitive

    nmlfi3_efold_primitive = nmlfi_efold_primitive(x,xi,p)

  end function nmlfi3_efold_primitive



  function nmlfi3_parametric_hbar_trajectory(bfold,hbarend,xi,p)
    implicit none
    real(kp), intent(in) :: bfold, hbarend, xi,p
    real(kp) :: nmlfi3_parametric_hbar_trajectory
    real(kp) :: hbarmin, hbarmax

    if (.not.nmlfi3_check_params(xi,p)) then
       write(*,*)'xi= p= ',xi,p
       stop 'nmlfi3_parametric_hbar_trajectory: NMLFI3 does not exist for these parameters!'
    endif
       
    hbarmin = nmlfi_hbar_potmax(xi,p)

    hbarmax = hbarBig
        
    nmlfi3_parametric_hbar_trajectory = nmlfi_parametric_hbar_trajectory(bfold,hbarend,hbarmin,hbarmax,xi,p)

  end function nmlfi3_parametric_hbar_trajectory

  
!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function nmlfi3_x_trajectory(bfold,xend,xi,p)
    implicit none
    real(kp), intent(in) :: bfold, xend, xi,p
    real(kp) :: nmlfi3_x_trajectory
    real(kp) :: hbar, hbarend

    hbarend = nmlfi_hbar(xend,xi)
    
    hbar = nmlfi3_parametric_hbar_trajectory(bfold,hbarend,xi,p)

    nmlfi3_x_trajectory = nmlfi_x(hbar,xi)
            
  end function nmlfi3_x_trajectory


end module nmlfi3sr
