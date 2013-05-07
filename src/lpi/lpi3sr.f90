!slow-roll functions for the logarithmic potential 3 inflation models
!
!V(phi) = M**4 x**p log(x)**q
!
!x = phi/phi0
!phi0 = phi0/Mp
!0 < x < xpotmax

module lpi3sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent


  use lpicommon, only : lpi_norm_potential, lpi_epsilon_one, lpi_epsilon_two
  use lpicommon, only : lpi_epsilon_three,  lpi_norm_deriv_potential 
  use lpicommon, only : lpi_norm_deriv_second_potential, lpi_x_potmax
  use lpicommon, only : lpi_efold_primitive, find_lpi_x_trajectory
  use lpicommon, only : lpi23_sanity_check

  implicit none

  private

 
  public lpi3_norm_potential, lpi3_epsilon_one, lpi3_epsilon_two, lpi3_epsilon_three
  public lpi3_norm_deriv_potential, lpi3_norm_deriv_second_potential
  public lpi3_x_endinf, lpi3_x_trajectory, lpi3_efold_primitive


 
contains



!returns V/M**4
  function lpi3_norm_potential(x,p,q,phi0)
    implicit none
    real(kp) :: lpi3_norm_potential
    real(kp), intent(in) :: x,p,q
    real(kp), intent(in) :: phi0

    call lpi23_sanity_check(q=q)

    lpi3_norm_potential = lpi_norm_potential(x,p,q)

  end function lpi3_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function lpi3_norm_deriv_potential(x,p,q,phi0)
    implicit none
    real(kp) :: lpi3_norm_deriv_potential
    real(kp), intent(in) :: x,p,q
    real(kp), intent(in) :: phi0

    call lpi23_sanity_check(q=q)

    lpi3_norm_deriv_potential = lpi_norm_deriv_potential(x,p,q)

  end function lpi3_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function lpi3_norm_deriv_second_potential(x,p,q,phi0)
    implicit none
    real(kp) :: lpi3_norm_deriv_second_potential
    real(kp), intent(in) :: x,p,q
    real(kp), intent(in) :: phi0

    call lpi23_sanity_check(q=q)

    lpi3_norm_deriv_second_potential = lpi_norm_deriv_second_potential(x,p,q)

  end function lpi3_norm_deriv_second_potential



!epsilon_one(x)
  function lpi3_epsilon_one(x,p,q,phi0)    
    implicit none
    real(kp) :: lpi3_epsilon_one
    real(kp), intent(in) :: x,p,q
    real(kp), intent(in) :: phi0

    call lpi23_sanity_check(q=q)

    lpi3_epsilon_one = lpi_epsilon_one(x,p,q,phi0)
    
  end function lpi3_epsilon_one


!epsilon_two(x)
  function lpi3_epsilon_two(x,p,q,phi0)    
    implicit none
    real(kp) :: lpi3_epsilon_two
    real(kp), intent(in) :: x,p,q
    real(kp), intent(in) :: phi0

    call lpi23_sanity_check(q=q)

    lpi3_epsilon_two = lpi_epsilon_two(x,p,q,phi0)
    
  end function lpi3_epsilon_two


!epsilon_three(x)
  function lpi3_epsilon_three(x,p,q,phi0)    
    implicit none
    real(kp) :: lpi3_epsilon_three
    real(kp), intent(in) :: x,p,q
    real(kp), intent(in) :: phi0

    call lpi23_sanity_check(q=q)
    
    lpi3_epsilon_three = lpi_epsilon_three(x,p,q,phi0)

    
  end function lpi3_epsilon_three



!returns the value for x=phi/Mp defined as epsilon1=1, where inflation ends
  function lpi3_x_endinf(p,q,phi0)
    implicit none
    real(kp), intent(in) :: phi0,p,q
    real(kp) :: lpi3_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lpi3Data

    call lpi23_sanity_check(q=q)

    mini = epsilon(1._kp)
    maxi = lpi_x_potmax(p,q)

    lpi3Data%real1 = p
    lpi3Data%real2 = q
    lpi3Data%real3 = phi0
    
    lpi3_x_endinf = zbrent(find_lpi3_x_endinf,mini,maxi,tolFind,lpi3Data)
      
  end function lpi3_x_endinf


 
  function find_lpi3_x_endinf(x,lpiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lpiData
    real(kp) :: find_lpi3_x_endinf
    real(kp) :: phi0,p,q

    p = lpiData%real1
    q = lpiData%real2
    phi0 = lpiData%real3

    find_lpi3_x_endinf = lpi3_epsilon_one(x,p,q,phi0)-1._kp
   
  end function find_lpi3_x_endinf



  !this is integral[V(phi)/V'(phi) dphi]
  function lpi3_efold_primitive(x,p,q,phi0)
    implicit none
    real(kp), intent(in) :: x,phi0,p,q
    real(kp) :: lpi3_efold_primitive

    call lpi23_sanity_check(q=q)

    lpi3_efold_primitive = lpi_efold_primitive(x,p,q,phi0)

  end function lpi3_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function lpi3_x_trajectory(bfold,xend,p,q,phi0)
    implicit none
    real(kp), intent(in) :: bfold, phi0,p,q, xend
    real(kp) :: lpi3_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,xpotmax
    type(transfert) :: lpi3Data

    call lpi23_sanity_check(q=q)

    xpotmax = lpi_x_potmax(p,q)

    if (xend.gt.xpotmax) then
       write(*,*)'xend= xpotmax= ',xend,xpotmax
       stop 'lpi3_x_trajectory: xend should be < xpotmax'
    endif
  
    mini = xend
    maxi = xpotmax

    lpi3Data%real1 = p
    lpi3Data%real2 = q
    lpi3Data%real3 = phi0
    lpi3Data%real4 = -bfold + lpi_efold_primitive(xend,p,q,phi0)
    
    lpi3_x_trajectory = zbrent(find_lpi_x_trajectory,mini,maxi,tolFind,lpi3Data)
       
  end function lpi3_x_trajectory


  
end module lpi3sr
