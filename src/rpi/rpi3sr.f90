!slow-roll functions for the R-R^p inflation potential
!in the rpi3 regime, i.e. for p<1
!
!V(phi) = M^4 * exp(-2y) * [ exp(y) - 1 ]^(2p/(2p-1))
!
!y = phi/Mp * sqrt(2/3)


module rpi3sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use rpicommon, only : rpi_norm_potential, rpi_norm_deriv_potential
  use rpicommon, only : rpi_norm_deriv_second_potential
  use rpicommon, only : rpi_epsilon_one, rpi_epsilon_two, rpi_epsilon_three
  use rpicommon, only : rpih_efold_primitive, rpih_x_trajectory
  use rpicommon, only : rpi_efold_primitive, find_rpi_x_trajectory
  use rpicommon, only : rpi_x_epsoneunity,rpiBig
  implicit none

  private

  public rpi3_norm_potential, rpi3_epsilon_one, rpi3_epsilon_two
  public rpi3_epsilon_three
  public rpi3_x_endinf, rpi3_efold_primitive, rpi3_x_trajectory
  public rpi3_norm_deriv_potential, rpi3_norm_deriv_second_potential 

 
contains
!returns V/M^4
  function rpi3_norm_potential(y,p)
    implicit none
    real(kp) :: rpi3_norm_potential
    real(kp), intent(in) :: y,p

    rpi3_norm_potential = rpi_norm_potential(y,p)

  end function rpi3_norm_potential



!returns the first derivative of the potential with respect to y,
!divided by M^4
  function rpi3_norm_deriv_potential(y,p)
    implicit none
    real(kp) :: rpi3_norm_deriv_potential
    real(kp), intent(in) :: y,p

    rpi3_norm_deriv_potential = rpi_norm_deriv_potential(y,p)

  end function rpi3_norm_deriv_potential



!returns the second derivative of the potential with respect to y,
!divided by M^4
  function rpi3_norm_deriv_second_potential(y,p)
    implicit none
    real(kp) :: rpi3_norm_deriv_second_potential
    real(kp), intent(in) :: y,p

    rpi3_norm_deriv_second_potential = rpi_norm_deriv_second_potential(y,p)

  end function rpi3_norm_deriv_second_potential



!epsilon_one(y)
  function rpi3_epsilon_one(y,p)    
    implicit none
    real(kp) :: rpi3_epsilon_one
    real(kp), intent(in) :: y,p


    rpi3_epsilon_one = rpi_epsilon_one(y,p)


  end function rpi3_epsilon_one


!epsilon_two(y)
  function rpi3_epsilon_two(y,p)    
    implicit none
    real(kp) :: rpi3_epsilon_two
    real(kp), intent(in) :: y,p

    rpi3_epsilon_two = rpi_epsilon_two(y,p)

  end function rpi3_epsilon_two


!epsilon_three(y)
  function rpi3_epsilon_three(y,p)    
    implicit none
    real(kp) :: rpi3_epsilon_three
    real(kp), intent(in) :: y,p

    rpi3_epsilon_three = rpi_epsilon_three(y,p)

  end function rpi3_epsilon_three


 
  function rpi3_x_endinf(p)
    implicit none
    real(kp), intent(in) :: p
    real(kp) :: rpi3_x_endinf

    rpi3_x_endinf = rpi_x_epsoneunity(p)
  
  end function rpi3_x_endinf


 
  !this is integral[V(phi)/V'(phi) dphi]
  function rpi3_efold_primitive(y,p)
    implicit none
    real(kp), intent(in) :: y,p
    real(kp) :: rpi3_efold_primitive

    if (p.eq.1._kp) then
       rpi3_efold_primitive = rpih_efold_primitive(y,p)
       return
    endif


    rpi3_efold_primitive = rpi_efold_primitive(y,p)

  end function rpi3_efold_primitive


  !returns y at bfold=-efolds before the end of inflation, ie N-Nend
  function rpi3_x_trajectory(bfold,yend,p)
    implicit none
    real(kp), intent(in) :: bfold, p, yend
    real(kp) :: rpi3_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rpi3Data


    if (p.eq.1._kp) then !Higgs Inflation Model (HI)
       rpi3_x_trajectory = rpih_x_trajectory(bfold,yend,p)
       return
    endif
    
    mini = yEnd
    maxi = rpiBig

    rpi3Data%real1 = p
    rpi3Data%real2 = -bfold + rpi3_efold_primitive(yend,p)

    rpi3_x_trajectory = zbrent(find_rpi_x_trajectory,mini,maxi,tolFind,rpi3Data)

  end function rpi3_x_trajectory


end module rpi3sr
