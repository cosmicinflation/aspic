!slow-roll functions for the R + R^2/m^2 + alpha R^3/m^4 inflation potential
!with alpha > 0 and x < xVmax
!
!V(m=k2) = 
!
!x = phi/Mp
!dx/dk = [4 sqrt(2)/pi] sqrt{K(k2) K(1-k2)}/k2

module disr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public di_norm_potential, di_epsilon_one, di_epsilon_two, di_epsilon_three
  public di_x_endinf, di_efold_primitive, di_x_trajectory
  public di_norm_deriv_potential, di_norm_deriv_second_potential 
  public di_check_params, di_numacc_xinimax

contains


  function di_check_params(f)
    implicit none
    logical :: di_check_params
    real(kp), intent(in) :: f

    di_check_params = (alpha.ge.0._kp)

  end function di_check_params



  function di_norm_potential(x,f)
    implicit none
    real(kp) :: di_norm_potential
    real(kp), intent(in) :: x,f

    real(kp) :: k2

    k = di_k(x)
    k2 = k*k
    
    di_norm_potential = di_norm_parametric_potential(k2)

  end function di_norm_potential




  function di_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: di_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

    di_norm_deriv_potential = ccsi_norm_deriv_potential(x,alpha)

  end function di_norm_deriv_potential




  function di_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: di_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    di_norm_deriv_second_potential = ccsi_norm_deriv_second_potential(x,alpha)

  end function di_norm_deriv_second_potential


  
  function di_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: di_epsilon_one
    real(kp), intent(in) :: x,alpha

    di_epsilon_one = ccsi_epsilon_one(x,alpha)

  end function di_epsilon_one



  function di_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: di_epsilon_two
    real(kp), intent(in) :: x,alpha

    di_epsilon_two = ccsi_epsilon_two(x,alpha)

  end function di_epsilon_two



  function di_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: di_epsilon_three
    real(kp), intent(in) :: x,alpha

    di_epsilon_three = ccsi_epsilon_three(x,alpha)

  end function di_epsilon_three


 
  function di_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: di_x_endinf

    real(kp), dimension(2) :: xepsone

    xepsone = ccsi_x_epsoneunity(alpha)

    di_x_endinf = xepsone(1)
  
  end function di_x_endinf

!return the maximal positive value of xini for ensuring eps1 >
!numerical accuracy, this is xVmax - smallterm
  function di_numacc_xinimax(alpha)
    implicit none
    real(kp) :: di_numacc_xinimax
    real(kp), intent(in) :: alpha

    real(kp), dimension(2) :: xnumacc

    if (.not.di_check_params(alpha)) then
       stop 'di_numacc_xinimax: di requires alpha=>0'
    endif
      
    xnumacc = ccsi_numacc_x_epsonenull(alpha)

    di_numacc_xinimax = xnumacc(1)
   
    
  end function di_numacc_xinimax
 
  !this is integral[V(phi)/V'(phi) dphi]
  function di_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: di_efold_primitive

    real(kp) :: xVmax
   
    if (alpha.eq.0._kp) then
       di_efold_primitive = ccsih_efold_primitive(x,alpha)
       return
    endif

    xVmax = ccsi_x_potmax(alpha)

    if (x.gt.xVmax) stop 'di_efold_primitive: x > xVmax'

    di_efold_primitive = ccsi_efold_primitive(x,alpha)

  end function di_efold_primitive


  !returns y at bfold=-efolds before the end of inflation, ie N-Nend
  function di_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: di_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: diData
    real(kp) :: xVmax

    if (.not.di_check_params(alpha)) then
       stop 'di_x_trajectory: di requires alpha=>0'
    endif

    if (alpha.eq.0._kp) then !Higgs Inflation Model (HI)
       di_x_trajectory = ccsih_x_trajectory(bfold,xend,alpha)
       return
    endif

    xVmax = ccsi_x_potmax(alpha)

    if (xend.gt.xVmax) stop 'di_x_trajectory: xend > xVmax'
   
    mini = xEnd
    maxi = di_numacc_xinimax(alpha)

    diData%real1 = alpha
    diData%real2 = -bfold + di_efold_primitive(xend,alpha)

    di_x_trajectory = zbrent(find_ccsi_x_trajectory,mini,maxi,tolFind,diData)

  end function di_x_trajectory


end module disr
