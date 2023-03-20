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
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use nmlficommon, only : hbarBig, pplus, nmlfi_hbar, nmlfi_x, nmlfi_xizero
  use nmlficommon, only : nmlfi_norm_potential, nmlfi_norm_deriv_potential, nmlfi_norm_deriv_second_potential
  use nmlficommon, only : nmlfi_epsilon_one, nmlfi_epsilon_two, nmlfi_epsilon_three
  use nmlficommon, only : nmlfi_efold_primitive, nmlfi_parametric_hbar_trajectory
  use nmlficommon, only : nmlfi_hbarsquare_epsoneunity, nmlfi_hbar_potmax
  use nmlficommon, only : nmlfi_numacc_hbarsquare_epsonenull, nmlfi_parametric_efold_primitive
  
  implicit none

  private

  public nmlfi1_norm_potential, nmlfi1_norm_deriv_potential, nmlfi1_norm_deriv_second_potential
  public nmlfi1_epsilon_one, nmlfi1_epsilon_two, nmlfi1_epsilon_three
  public nmlfi1_x_endinf, nmlfi1_efold_primitive, nmlfi1_x_trajectory, nmlfi1_ximax
  public nmlfi1_hbar_endinf, nmlfi1_parametric_hbar_trajectory, nmlfi1_check_params

!for p<4 when nmlfi1 is confined at x < xVmax, some numerical issues
!occur at the top of the potential
  public nmlfi1_numacc_hbarinimax, nmlfi1_numacc_xinimax, nmlfi1_numacc_efoldmax
  public nmlfi1_numacc_ximax
  

  
contains


  function nmlfi1_check_params(xi,p)
    implicit none
    logical :: nmlfi1_check_params
    real(kp), intent(in) :: xi,p

    nmlfi1_check_params = (p.lt.pplus) &
         .or. ( (p.ge.pplus).and.(xi.lt.nmlfi_xizero(p)) )

  end function nmlfi1_check_params


  
  function nmlfi1_norm_potential(x,xi,p)
    implicit none
    real(kp) :: nmlfi1_norm_potential
    real(kp), intent(in) :: x,xi,p
    
    nmlfi1_norm_potential = nmlfi_norm_potential(x,xi,p)

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

  
  function nmlfi1_epsilon_one(x,xi,p)
    implicit none
    real(kp) :: nmlfi1_epsilon_one
    real(kp), intent(in) :: x,xi,p

    nmlfi1_epsilon_one = nmlfi_epsilon_one(x,xi,p)

  end function nmlfi1_epsilon_one


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
    real(kp), intent(in) :: xi,p
    
    real(kp) :: hbarendinf
    
    hbarendinf = nmlfi1_hbar_endinf(xi,p)

    nmlfi1_x_endinf = nmlfi_x(hbarendinf,xi)
    
  end function nmlfi1_x_endinf


  
  function nmlfi1_hbar_endinf(xi,p)
    implicit none
    real(kp) :: nmlfi1_hbar_endinf
    real(kp), intent(in) :: xi,p

    real(kp), dimension(2) :: hbarepsone2

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


!no inflation for p > pplus and xi > xizero  
  function nmlfi1_ximax(p)
    implicit none
    real(kp) :: nmlfi1_ximax
    real(kp), intent(in) :: p

    if (p.ge.pplus) then
       nmlfi1_ximax = nmlfi_xizero(p)
    else
       nmlfi1_ximax = huge(1._kp)
    endif

  end function nmlfi1_ximax
  


!for p<4, an infinite number of e-folds can be made close to the top
!of the potential, but this triggers numerical issues as eps1->0, plus
!ns->1 and this is not of immediate interest. This function returns
!the maximal possible value of hbarini, close to the top of the
!potential, starting at eps1=machine precision
  function nmlfi1_numacc_hbarinimax(xi,p)
    implicit none
    real(kp) :: nmlfi1_numacc_hbarinimax
    real(kp), intent(in) :: xi,p

    real(kp), dimension(2) :: hbar2null

!the potential has no maximum, we return a large number    
    if (p.ge.4._kp) then
       nmlfi1_numacc_hbarinimax = hbarBig
       return
    endif
    
    hbar2null = nmlfi_numacc_hbarsquare_epsonenull(xi,p)

    if ((hbar2null(1).le.0._kp).and.(p.lt.4._kp)) then
       write(*,*)'xi= p= hbar2null= ',xi,p,hbar2null(:)
       stop ' nmlfi1_numacc_hbarinimax: bug!?'
    endif

    nmlfi1_numacc_hbarinimax = sqrt(hbar2null(1))

!safety check
    if (nmlfi1_numacc_hbarinimax.gt.nmlfi_hbar_potmax(xi,p)) then
       stop 'nmlfi1_numacc_hbarinimax: hbarinimax > hbarpotmax!?'
    endif

  end function nmlfi1_numacc_hbarinimax
    
  

  function nmlfi1_numacc_xinimax(xi,p)
    implicit none
    real(kp) :: nmlfi1_numacc_xinimax
    real(kp), intent(in) :: xi,p
    real(kp) :: hbarinimax

    hbarinimax = nmlfi1_numacc_hbarinimax(xi,p)

    nmlfi1_numacc_xinimax = nmlfi_x(hbarinimax,xi)
    
  end function nmlfi1_numacc_xinimax


!the maximal number of efold computable at the current numerical accuracy
  function nmlfi1_numacc_efoldmax(xi,p)
    implicit none
    real(kp) :: nmlfi1_numacc_efoldmax
    real(kp), intent(in) :: xi,p

    real(kp) :: hbarend, hbarinimax
    
    hbarend = nmlfi1_hbar_endinf(xi,p)
    hbarinimax = nmlfi1_numacc_hbarinimax(xi,p)

    nmlfi1_numacc_efoldmax = -nmlfi_parametric_efold_primitive(hbarend,xi,p) &
         + nmlfi_parametric_efold_primitive(hbarinimax,xi,p)

  end function nmlfi1_numacc_efoldmax


!the number of-efolds at a dependence in 1/xi such that, at fixed
!numerical accuracy, getting xi small gives more e-folds. This
!function returns the maximal value of xi to get efold
  function nmlfi1_numacc_ximax(efold,p)
    implicit none
    real(kp) :: nmlfi1_numacc_ximax
    real(kp), intent(in) :: efold,p

    type(transfert) :: nmlfi1Data
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi

    real(kp), save :: efoldSave = 0._kp
    real(kp), save :: pSave = 0._kp
    real(kp), save :: ximaxSave = 0._kp
!$omp threadprivate(efoldSave,ximaxSave)

    
!accelerator for repeated calls
    if ((efold.eq.efoldSave).and.(p.eq.pSave)) then
       nmlfi1_numacc_ximax = ximaxSave
       return
    else
       efoldSave = efold
       pSave = p
    endif        

    mini = epsilon(1._kp)
    maxi = nmlfi1_ximax(p) - epsilon(1._kp)

    nmlfi1Data%real1 = efold
    nmlfi1Data%real2 = p
    
    nmlfi1_numacc_ximax = zbrent(find_nmlfi1_ximax,mini,maxi,tolFind,nmlfi1Data)

    ximaxSave = nmlfi1_numacc_ximax
    
  end function nmlfi1_numacc_ximax



  function find_nmlfi1_ximax(xi,nmlfi1Data)
    implicit none
    real(kp), intent(in) :: xi
    type(transfert), optional, intent(inout) :: nmlfi1Data
    real(kp) :: find_nmlfi1_ximax
    real(kp) :: efold,p

    efold = nmlfi1Data%real1
    p = nmlfi1Data%real2
    
    find_nmlfi1_ximax = efold  - nmlfi1_numacc_efoldmax(xi,p)

  end function find_nmlfi1_ximax

  
  
  
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
    real(kp) :: hbar, hbarend

    hbarend = nmlfi_hbar(xend,xi)
    
    hbar = nmlfi1_parametric_hbar_trajectory(bfold,hbarend,xi,p)

    nmlfi1_x_trajectory = nmlfi_x(hbar,xi)
            
  end function nmlfi1_x_trajectory


end module nmlfi1sr
