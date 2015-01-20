!slow-roll functions for the R + R^2/m^2 + alpha R^3/m^4 inflation potential
!with alpha < 0
!
!V(phi) = M^4 * exp(-2x) * [ exp(x) - 1 ]^2/{1 + sqrt[1+3 alpha(exp(x)-1)]}^3 
!  * {1 + sqrt[1 + 3 alpha (exp(x)-1)] + 2 alpha (exp(x)-1) }
!
!x = phi/Mp * sqrt(2/3)


module ccsi3sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ccsicommon, only : ccsi_norm_potential, ccsi_norm_deriv_potential
  use ccsicommon, only : ccsi_norm_deriv_second_potential
  use ccsicommon, only : ccsi_epsilon_one, ccsi_epsilon_two, ccsi_epsilon_three
  use ccsicommon, only : ccsih_efold_primitive, ccsih_x_trajectory, ccsi_x_potmax
  use ccsicommon, only : ccsi_efold_primitive, find_ccsi_x_trajectory
  use ccsicommon, only : ccsi_x_epsoneunity, ccsi_alphamin, ccsi_xmax
  use ccsicommon, only : ccsiBig
  implicit none

  private

  public ccsi3_norm_potential, ccsi3_epsilon_one, ccsi3_epsilon_two, ccsi3_epsilon_three
  public ccsi3_efold_primitive, ccsi3_x_trajectory, ccsi3_x_endinf
  public ccsi3_norm_deriv_potential, ccsi3_norm_deriv_second_potential 
  public ccsi3_check_params, ccsi3_xinimax, ccsi3_alphamin

contains


  function ccsi3_check_params(alpha)
    implicit none
    logical :: ccsi3_check_params
    real(kp), intent(in) :: alpha

    ccsi3_check_params = (alpha.lt.0._kp).and.(alpha.gt.ccsi_alphamin())

  end function ccsi3_check_params



  function ccsi3_norm_potential(x,alpha)
    implicit none
    real(kp) :: ccsi3_norm_potential
    real(kp), intent(in) :: x,alpha

    ccsi3_norm_potential = ccsi_norm_potential(x,alpha)

  end function ccsi3_norm_potential




  function ccsi3_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: ccsi3_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

    ccsi3_norm_deriv_potential = ccsi_norm_deriv_potential(x,alpha)

  end function ccsi3_norm_deriv_potential




  function ccsi3_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: ccsi3_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    ccsi3_norm_deriv_second_potential = ccsi_norm_deriv_second_potential(x,alpha)

  end function ccsi3_norm_deriv_second_potential


  
  function ccsi3_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: ccsi3_epsilon_one
    real(kp), intent(in) :: x,alpha

    ccsi3_epsilon_one = ccsi_epsilon_one(x,alpha)

  end function ccsi3_epsilon_one



  function ccsi3_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: ccsi3_epsilon_two
    real(kp), intent(in) :: x,alpha

    ccsi3_epsilon_two = ccsi_epsilon_two(x,alpha)

  end function ccsi3_epsilon_two



  function ccsi3_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: ccsi3_epsilon_three
    real(kp), intent(in) :: x,alpha

    ccsi3_epsilon_three = ccsi_epsilon_three(x,alpha)

  end function ccsi3_epsilon_three


  function ccsi3_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: ccsi3_x_endinf

    real(kp), dimension(2) :: xepsone

    xepsone = ccsi_x_epsoneunity(alpha)

    ccsi3_x_endinf = xepsone(1)
  
  end function ccsi3_x_endinf


!return the maximal positive value of xini such that the potential is
!still defined (xini < xmax) and epsilon1 < 1
  function ccsi3_xinimax(alpha)
    implicit none
    real(kp) :: ccsi3_xinimax
    real(kp), intent(in) :: alpha
    real(kp) :: xmax
    real(kp), dimension(2) :: xepsone

    if (.not.ccsi3_check_params(alpha)) then
       stop 'ccsi3_xinimax: ccsi3 requires alpha<0'
    endif

    xmax = ccsi_xmax(alpha)
    xepsone = ccsi_x_epsoneunity(alpha)

!this cannot happen
    if (xepsone(2).gt.xmax) stop 'ccsi3_xinimax: internal error'

    ccsi3_xinimax = xepsone(2)
       
  end function ccsi3_xinimax
 

!returns the minimum (<0) value of alpha to get at most efoldMax
!of inflation (from xinimax to xend)
  function ccsi3_alphamin(efoldMax)
    implicit none
    real(kp) :: ccsi3_alphamin
    real(kp), intent(in) :: efoldMax
    type(transfert) :: ccsi3Data
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini, maxi

    mini = ccsi_alphamin() + epsilon(1._kp)
    maxi = -epsilon(1._kp)

    ccsi3Data%real1 = efoldMax
    ccsi3_alphamin = zbrent(find_ccsi3_alphamin,mini,maxi,tolFind,ccsi3Data)

  end function ccsi3_alphamin

  
  function find_ccsi3_alphamin(alpha,ccsi3Data)
    implicit none
    real(kp), intent(in) :: alpha
    type(transfert), optional, intent(inout) :: ccsi3Data
    real(kp) :: find_ccsi3_alphamin
    real(kp) :: xinimax, xend, efoldMax

    efoldMax = ccsi3Data%real1
    xinimax = ccsi3_xinimax(alpha)
    xend = ccsi3_x_endinf(alpha)

    find_ccsi3_alphamin = efoldMax + ccsi3_efold_primitive(xend,alpha) &
         - ccsi3_efold_primitive(xinimax,alpha)

  end function find_ccsi3_alphamin



!this is integral[V(phi)/V'(phi) dphi]
  function ccsi3_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: ccsi3_efold_primitive
   
    if (alpha.eq.0._kp) then
       ccsi3_efold_primitive = ccsih_efold_primitive(x,alpha)
       return
    endif

    ccsi3_efold_primitive = ccsi_efold_primitive(x,alpha)

  end function ccsi3_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ccsi3_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: ccsi3_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi, efoldMax, xinimax
    type(transfert) :: ccsi3Data

    if (.not.ccsi3_check_params(alpha)) then
       stop 'ccsi3_x_trajectory: ccsi3 requires alpha<0'
    endif

    if (alpha.eq.0._kp) then !Higgs Inflation Model (HI)
       ccsi3_x_trajectory = ccsih_x_trajectory(bfold,xend,alpha)
       return
    endif

    xinimax = ccsi3_xinimax(alpha)
            
    efoldMax = -ccsi3_efold_primitive(xend,alpha) &
         + ccsi3_efold_primitive(xinimax,alpha)

    if (-bfold.gt.efoldMax) then
       write(*,*)'ccsi3_x_trajectory: not enough efolds!'
       write(*,*)'efold requested= efold maxi= ',-bfold,efoldMax       
       stop
    endif

    mini = xEnd
    maxi = xinimax
       
    ccsi3Data%real1 = alpha
    ccsi3Data%real2 = -bfold + ccsi3_efold_primitive(xend,alpha)

    ccsi3_x_trajectory = zbrent(find_ccsi_x_trajectory,mini,maxi,tolFind,ccsi3Data)

  end function ccsi3_x_trajectory


end module ccsi3sr
