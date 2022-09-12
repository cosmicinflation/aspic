!slow-roll functions for the R + R^2/m^2 + alpha R^3/m^4 inflation potential
!with alpha > 0 and x < xVmax
!
!V(phi) = M^4 * exp(-2x) * [ exp(x) - 1 ]^2/{1 + sqrt[1+3 alpha(exp(x)-1)]}^3 
!  * {1 + sqrt[1 + 3 alpha (exp(x)-1)] + 2 alpha (exp(x)-1) }
!
!x = phi/Mp * sqrt(2/3)


module ccsi1sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ccsicommon, only : ccsi_norm_potential, ccsi_norm_deriv_potential
  use ccsicommon, only : ccsi_norm_deriv_second_potential
  use ccsicommon, only : ccsi_epsilon_one, ccsi_epsilon_two, ccsi_epsilon_three
  use ccsicommon, only : ccsih_efold_primitive, ccsih_x_trajectory, ccsi_x_potmax
  use ccsicommon, only : ccsi_efold_primitive, find_ccsi_x_trajectory
  use ccsicommon, only : ccsi_x_epsoneunity, ccsi_numacc_x_epsonenull
  implicit none

  private

  public ccsi1_norm_potential, ccsi1_epsilon_one, ccsi1_epsilon_two, ccsi1_epsilon_three
  public ccsi1_x_endinf, ccsi1_efold_primitive, ccsi1_x_trajectory
  public ccsi1_norm_deriv_potential, ccsi1_norm_deriv_second_potential 
  public ccsi1_check_params, ccsi1_numacc_xinimax, ccsi1_numacc_efoldmax
  public ccsi1_numacc_alphamax
contains


  function ccsi1_check_params(alpha)
    implicit none
    logical :: ccsi1_check_params
    real(kp), intent(in) :: alpha

    ccsi1_check_params = (alpha.ge.0._kp)

  end function ccsi1_check_params



  function ccsi1_norm_potential(x,alpha)
    implicit none
    real(kp) :: ccsi1_norm_potential
    real(kp), intent(in) :: x,alpha

    ccsi1_norm_potential = ccsi_norm_potential(x,alpha)

  end function ccsi1_norm_potential




  function ccsi1_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: ccsi1_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

    ccsi1_norm_deriv_potential = ccsi_norm_deriv_potential(x,alpha)

  end function ccsi1_norm_deriv_potential




  function ccsi1_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: ccsi1_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    ccsi1_norm_deriv_second_potential = ccsi_norm_deriv_second_potential(x,alpha)

  end function ccsi1_norm_deriv_second_potential


  
  function ccsi1_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: ccsi1_epsilon_one
    real(kp), intent(in) :: x,alpha

    ccsi1_epsilon_one = ccsi_epsilon_one(x,alpha)

  end function ccsi1_epsilon_one



  function ccsi1_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: ccsi1_epsilon_two
    real(kp), intent(in) :: x,alpha

    ccsi1_epsilon_two = ccsi_epsilon_two(x,alpha)

  end function ccsi1_epsilon_two



  function ccsi1_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: ccsi1_epsilon_three
    real(kp), intent(in) :: x,alpha

    ccsi1_epsilon_three = ccsi_epsilon_three(x,alpha)

  end function ccsi1_epsilon_three


 
  function ccsi1_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: ccsi1_x_endinf

    real(kp), dimension(2) :: xepsone

    xepsone = ccsi_x_epsoneunity(alpha)

    ccsi1_x_endinf = xepsone(1)
  
  end function ccsi1_x_endinf

!return the maximal positive value of xini for ensuring eps1 >
!numerical accuracy, this is xVmax - smallterm
  function ccsi1_numacc_xinimax(alpha)
    implicit none
    real(kp) :: ccsi1_numacc_xinimax
    real(kp), intent(in) :: alpha

    real(kp), dimension(2) :: xnumacc

    if (.not.ccsi1_check_params(alpha)) then
       stop 'ccsi1_numacc_xinimax: ccsi1 requires alpha=>0'
    endif
      
    xnumacc = ccsi_numacc_x_epsonenull(alpha)

    ccsi1_numacc_xinimax = xnumacc(1)
   
    
  end function ccsi1_numacc_xinimax




!returns the minimum (<0) value of alpha to get at most efoldMax
!of inflation (from xinimax to xend)
  function ccsi1_numacc_alphamax(efoldMax)
    implicit none
    real(kp) :: ccsi1_numacc_alphamax
    real(kp), intent(in) :: efoldMax
    type(transfert) :: ccsi1Data
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini, maxi

    real(kp), save :: efoldSave = 0._kp
    real(kp), save :: alphamaxSave = 0._kp
!$omp threadprivate(efoldSave,alphamaxSave)
    
    if (efoldMax.eq.efoldSave) then
       ccsi1_numacc_alphamax = alphamaxSave
       return
    else
       efoldSave = efoldMax
    endif
    
    mini = epsilon(1._kp)
    maxi = 1._kp

    ccsi1Data%real1 = efoldMax
    ccsi1_numacc_alphamax  = zbrent(find_ccsi1_alphamax,mini,maxi,tolFind,ccsi1Data)

    alphamaxSave = ccsi1_numacc_alphamax
    
  end function ccsi1_numacc_alphamax

  
  function find_ccsi1_alphamax(alpha,ccsi1Data)
    implicit none
    real(kp), intent(in) :: alpha
    type(transfert), optional, intent(inout) :: ccsi1Data
    real(kp) :: find_ccsi1_alphamax
    real(kp) :: xinimax, xend, efoldMax

    efoldMax = ccsi1Data%real1

    find_ccsi1_alphamax = efoldMax  - ccsi1_numacc_efoldmax(alpha)

  end function find_ccsi1_alphamax
  
  

!return the maximal number of efold computable at this numerical accuracy
  function ccsi1_numacc_efoldmax(alpha)
    implicit none
    real(kp) :: ccsi1_numacc_efoldmax
    real(kp), intent(in) :: alpha

    real(kp) :: xendinf, xinimax
    
    xinimax = ccsi1_numacc_xinimax(alpha)
    xendinf = ccsi1_x_endinf(alpha)

    ccsi1_numacc_efoldmax = -ccsi1_efold_primitive(xendinf,alpha) &
         + ccsi1_efold_primitive(xinimax,alpha)

  end function ccsi1_numacc_efoldmax

  
  !this is integral[V(phi)/V'(phi) dphi]
  function ccsi1_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: ccsi1_efold_primitive

    real(kp) :: xVmax
   
    if (alpha.eq.0._kp) then
       ccsi1_efold_primitive = ccsih_efold_primitive(x,alpha)
       return
    endif

    xVmax = ccsi_x_potmax(alpha)

    if (x.gt.xVmax) stop 'ccsi1_efold_primitive: x > xVmax'

    ccsi1_efold_primitive = ccsi_efold_primitive(x,alpha)

  end function ccsi1_efold_primitive


  !returns y at bfold=-efolds before the end of inflation, ie N-Nend
  function ccsi1_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: ccsi1_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ccsi1Data
    real(kp) :: xVmax

    if (.not.ccsi1_check_params(alpha)) then
       stop 'ccsi1_x_trajectory: ccsi1 requires alpha=>0'
    endif

    if (alpha.eq.0._kp) then !Higgs Inflation Model (HI)
       ccsi1_x_trajectory = ccsih_x_trajectory(bfold,xend,alpha)
       return
    endif

    xVmax = ccsi_x_potmax(alpha)

    if (xend.gt.xVmax) stop 'ccsi1_x_trajectory: xend > xVmax'
   
    mini = xEnd
    maxi = ccsi1_numacc_xinimax(alpha)

    ccsi1Data%real1 = alpha
    ccsi1Data%real2 = -bfold + ccsi1_efold_primitive(xend,alpha)

    ccsi1_x_trajectory = zbrent(find_ccsi_x_trajectory,mini,maxi,tolFind,ccsi1Data)

  end function ccsi1_x_trajectory


end module ccsi1sr
