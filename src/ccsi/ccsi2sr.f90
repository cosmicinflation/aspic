!slow-roll functions for the R + R^2/m^2 + alpha R^3/m^4 inflation potential
!with alpha > 0 and x > xVmax
!
!V(phi) = M^4 * exp(-2x) * [ exp(x) - 1 ]^2/{1 + sqrt[1+3 alpha(exp(x)-1)]}^3 
!  * {1 + sqrt[1 + 3 alpha (exp(x)-1)] + 2 alpha (exp(x)-1) }
!
!x = phi/Mp * sqrt(2/3)


module ccsi2sr
  use infprec, only : kp,tolkp,toldp,transfert
  use inftools, only : zbrent
  use ccsicommon, only : ccsi_norm_potential, ccsi_norm_deriv_potential
  use ccsicommon, only : ccsi_norm_deriv_second_potential
  use ccsicommon, only : ccsi_epsilon_one, ccsi_epsilon_two, ccsi_epsilon_three
  use ccsicommon, only : ccsih_efold_primitive, ccsih_x_trajectory, ccsi_x_potmax
  use ccsicommon, only : ccsi_efold_primitive, find_ccsi_x_trajectory
  use ccsicommon, only : ccsi_x_epsoneunity, ccsi_numacc_x_epsonenull
  use ccsicommon, only : ccsiBig
  implicit none

  private

  public ccsi2_norm_potential, ccsi2_epsilon_one, ccsi2_epsilon_two, ccsi2_epsilon_three
  public ccsi2_efold_primitive, ccsi2_x_trajectory
  public ccsi2_norm_deriv_potential, ccsi2_norm_deriv_second_potential 
  public ccsi2_check_params, ccsi2_numacc_xinimin, ccsi2_numacc_xendmin

contains


  function ccsi2_check_params(alpha)
    implicit none
    logical :: ccsi2_check_params
    real(kp), intent(in) :: alpha

    ccsi2_check_params = (alpha.ge.0._kp)

  end function ccsi2_check_params



  function ccsi2_norm_potential(x,alpha)
    implicit none
    real(kp) :: ccsi2_norm_potential
    real(kp), intent(in) :: x,alpha

    ccsi2_norm_potential = ccsi_norm_potential(x,alpha)

  end function ccsi2_norm_potential




  function ccsi2_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: ccsi2_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

    ccsi2_norm_deriv_potential = ccsi_norm_deriv_potential(x,alpha)

  end function ccsi2_norm_deriv_potential




  function ccsi2_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: ccsi2_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    ccsi2_norm_deriv_second_potential = ccsi_norm_deriv_second_potential(x,alpha)

  end function ccsi2_norm_deriv_second_potential


  
  function ccsi2_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: ccsi2_epsilon_one
    real(kp), intent(in) :: x,alpha

    ccsi2_epsilon_one = ccsi_epsilon_one(x,alpha)

  end function ccsi2_epsilon_one



  function ccsi2_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: ccsi2_epsilon_two
    real(kp), intent(in) :: x,alpha

    ccsi2_epsilon_two = ccsi_epsilon_two(x,alpha)

  end function ccsi2_epsilon_two



  function ccsi2_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: ccsi2_epsilon_three
    real(kp), intent(in) :: x,alpha

    ccsi2_epsilon_three = ccsi_epsilon_three(x,alpha)

  end function ccsi2_epsilon_three


  

!return the minimal positive value of xini for ensuring eps1 >
!numerical accuracy, this is xVmax + smallterm
  function ccsi2_numacc_xinimin(alpha)
    implicit none
    real(kp) :: ccsi2_numacc_xinimin
    real(kp), intent(in) :: alpha

    real(kp), dimension(2) :: xnumacc

    if (.not.ccsi2_check_params(alpha)) then
       stop 'ccsi2_numacc_xinimax: ccsi2 requires alpha=>0'
    endif
      
    xnumacc = ccsi_numacc_x_epsonenull(alpha)

    ccsi2_numacc_xinimin = xnumacc(2)
   
  end function ccsi2_numacc_xinimin
 

!returns the minimal value of xend to get efold number of
!inflation while still starting in the domain xini > numacc_xinimin
   function ccsi2_numacc_xendmin(efold,alpha)
    implicit none
    real(kp) :: ccsi2_numacc_xendmin
    real(kp), intent(in) :: efold,alpha
    logical, parameter :: display = .false.
    real(kp) :: xinimin

    if (.not.ccsi2_check_params(alpha)) then
       stop 'ccsi2_numacc_xendmax: ccsi2 requires alpha>=0'
    endif

    xinimin = ccsi2_numacc_xinimin(alpha)

    ccsi2_numacc_xendmin = ccsi2_x_trajectory(efold,xinimin,alpha)

  end function ccsi2_numacc_xendmin



  !this is integral[V(phi)/V'(phi) dphi]
  function ccsi2_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: ccsi2_efold_primitive

    real(kp) :: xVmax
   
    if (alpha.eq.0._kp) then
       ccsi2_efold_primitive = ccsih_efold_primitive(x,alpha)
       return
    endif

    xVmax = ccsi_x_potmax(alpha)

    if (x.lt.xVmax) stop 'ccsi2_efold_primitive: x < xVmax'

    ccsi2_efold_primitive = ccsi_efold_primitive(x,alpha)

  end function ccsi2_efold_primitive


  !returns y at bfold=-efolds before the end of inflation, ie N-Nend
  function ccsi2_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: ccsi2_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ccsi2Data
    real(kp) :: xVmax

    if (.not.ccsi2_check_params(alpha)) then
       stop 'ccsi2_x_trajectory: ccsi2 requires alpha=>0'
    endif

    if (alpha.eq.0._kp) then !Higgs Inflation Model (HI)
       ccsi2_x_trajectory = ccsih_x_trajectory(bfold,xend,alpha)
       return
    endif

    xVmax = ccsi_x_potmax(alpha)

    if (xend.lt.xVmax) stop 'ccsi2_x_trajectory: xend < xVmax'
    
    if (bfold.le.0._kp) then
       mini = ccsi2_numacc_xinimin(alpha)
       maxi = xEnd
    else
       mini = xEnd
       maxi = ccsiBig
    endif


    ccsi2Data%real1 = alpha
    ccsi2Data%real2 = -bfold + ccsi2_efold_primitive(xend,alpha)

    ccsi2_x_trajectory = zbrent(find_ccsi_x_trajectory,mini,maxi,tolFind,ccsi2Data)

  end function ccsi2_x_trajectory


end module ccsi2sr
