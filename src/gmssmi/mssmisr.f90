!slow-roll functions for the MSSM inflation potential
!
!V(phi) = M^4 [ x^2 - 2/3 x^6 + 1/5 x^10 ]
!
!x = phi/Mp

module mssmisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
! this exists in recent fortran
! use specialinf, only : atanh
  use gmssmicommon, only : gmssmi_norm_potential, gmssmi_norm_deriv_potential
  use gmssmicommon, only : gmssmi_norm_deriv_second_potential
  use gmssmicommon, only : gmssmi_epsilon_one, gmssmi_epsilon_two, gmssmi_epsilon_three
  use gmssmicommon, only : gmssmi_x_endinf, gmssmi_x_epsonemin
  implicit none

  private

  public  mssmi_norm_potential, mssmi_epsilon_one, mssmi_epsilon_two, mssmi_epsilon_three
  public  mssmi_x_endinf, mssmi_efold_primitive, mssmi_x_trajectory
  public  mssmi_norm_deriv_potential, mssmi_norm_deriv_second_potential
  public  mssmi_x_epsonemin

 
contains


  function mssmi_alpha() !Returns alpha=1 which corresponds to the particular case MSSM inflation 
    implicit none
    real(kp) :: mssmi_alpha

    mssmi_alpha = 1._kp

  end function mssmi_alpha


!returns V/M**4
  function mssmi_norm_potential(x,phi0)
    implicit none
    real(kp) :: mssmi_norm_potential
    real(kp), intent(in) :: x,phi0

    mssmi_norm_potential = gmssmi_norm_potential(x,mssmi_alpha(),phi0)

  end function mssmi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function mssmi_norm_deriv_potential(x,phi0)
    implicit none
    real(kp) :: mssmi_norm_deriv_potential
    real(kp), intent(in) :: x,phi0

   mssmi_norm_deriv_potential = gmssmi_norm_deriv_potential(x,mssmi_alpha(),phi0)

  end function mssmi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function mssmi_norm_deriv_second_potential(x,phi0)
    implicit none
    real(kp) :: mssmi_norm_deriv_second_potential
    real(kp), intent(in) :: x,phi0

    mssmi_norm_deriv_second_potential = gmssmi_norm_deriv_second_potential(x,mssmi_alpha(),phi0)

  end function mssmi_norm_deriv_second_potential



!epsilon_one(x)
  function mssmi_epsilon_one(x,phi0)    
    implicit none
    real(kp) :: mssmi_epsilon_one
    real(kp), intent(in) :: x,phi0
    
    mssmi_epsilon_one = gmssmi_epsilon_one(x,mssmi_alpha(),phi0)
    
  end function mssmi_epsilon_one


!epsilon_two(x)
  function mssmi_epsilon_two(x,phi0)    
    implicit none
    real(kp) :: mssmi_epsilon_two
    real(kp), intent(in) :: x,phi0
    
    mssmi_epsilon_two = gmssmi_epsilon_two(x,mssmi_alpha(),phi0) 
    
  end function mssmi_epsilon_two


!epsilon_three(x)
  function mssmi_epsilon_three(x,phi0)    
    implicit none
    real(kp) :: mssmi_epsilon_three
    real(kp), intent(in) :: x,phi0
    
    mssmi_epsilon_three = gmssmi_epsilon_three(x,mssmi_alpha(),phi0) 
    
  end function mssmi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function mssmi_x_endinf(phi0)
    implicit none
    real(kp), intent(in) :: phi0
    real(kp) :: mssmi_x_endinf
    
    mssmi_x_endinf = gmssmi_x_endinf(mssmi_alpha(),phi0)
   

  end function mssmi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function mssmi_efold_primitive(x,phi0)
    implicit none
    real(kp), intent(in) :: x,phi0
    real(kp) :: mssmi_efold_primitive

    mssmi_efold_primitive = phi0**2*(x**2/20._kp-1._kp/15._kp*x**2/(x**4-1._kp) &
               + 2._kp/15._kp * atanh(x**2))

  end function mssmi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function mssmi_x_trajectory(bfold,xend,phi0)
    implicit none
    real(kp), intent(in) :: bfold, phi0, xend
    real(kp) :: mssmi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: mssmiData

  
    mini = xend
    maxi = mssmi_x_epsonemin()*(1._kp-epsilon(1._kp))

    mssmiData%real1 = phi0
    mssmiData%real2 = -bfold + mssmi_efold_primitive(xend,phi0)
    
    mssmi_x_trajectory = zbrent(find_mssmi_x_trajectory,mini,maxi,tolFind,mssmiData)
       
  end function mssmi_x_trajectory

  function find_mssmi_x_trajectory(x,mssmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: mssmiData
    real(kp) :: find_mssmi_x_trajectory
    real(kp) :: phi0,NplusNuend

    phi0= mssmiData%real1
    NplusNuend = mssmiData%real2

    find_mssmi_x_trajectory = mssmi_efold_primitive(x,phi0) - NplusNuend
   
  end function find_mssmi_x_trajectory



!Returns the position of the first local minimum of epsilon1
  function mssmi_x_epsonemin()   
    implicit none
    real(kp) :: mssmi_x_epsonemin

	mssmi_x_epsonemin = gmssmi_x_epsonemin(mssmi_alpha()) 
    
  end function mssmi_x_epsonemin

  
end module mssmisr
