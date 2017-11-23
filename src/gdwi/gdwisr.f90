!slow-roll functions for the generalized double well inflation potential
!
!V(phi) = M^4 * [(phi/phi0)^(2p) - 1]^2
!
!x=phi/phi0
!phi0=phi0/Mp


module gdwisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public  gdwi_norm_potential, gdwi_epsilon_one, gdwi_epsilon_two, gdwi_epsilon_three
  public  gdwi_x_endinf, gdwi_efold_primitive, gdwi_x_trajectory
  public  gdwi_norm_deriv_potential, gdwi_norm_deriv_second_potential
 
contains
!returns V/M^4
  function gdwi_norm_potential(x,p,phi0)
    implicit none
    real(kp) :: gdwi_norm_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: p,phi0

    gdwi_norm_potential = ((x*x)**p-1._kp)**2

  end function gdwi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function gdwi_norm_deriv_potential(x,p,phi0)
    implicit none
    real(kp) :: gdwi_norm_deriv_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: p,phi0

    gdwi_norm_deriv_potential = 4*p*x*(x*x)**(p-1)*((x*x)**p-1)
    
  end function gdwi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function gdwi_norm_deriv_second_potential(x,p,phi0)
    implicit none
    real(kp) :: gdwi_norm_deriv_second_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: p,phi0

    gdwi_norm_deriv_second_potential =  8*p**2*x*x*(x*x)**(-2 + 2*p) &
         + 8*(-1 + p)*p*x*x*(x*x)**(-2 + p)*(-1 + (x*x)**p) &
         + 4*p*(x*x)**(-1 + p)*(-1 + (x*x)**p)

  end function gdwi_norm_deriv_second_potential



!epsilon_one(x)
  function gdwi_epsilon_one(x,p,phi0)    
    implicit none
    real(kp) :: gdwi_epsilon_one
    real(kp), intent(in) :: x,p,phi0

  
    gdwi_epsilon_one = (8*p*p*(x*x)**(2*p-1))/(phi0**2*(-1+(x*x)**p)**2)

    
  end function gdwi_epsilon_one


!epsilon_two(x)
  function gdwi_epsilon_two(x,p,phi0)    
    implicit none
    real(kp) :: gdwi_epsilon_two
    real(kp), intent(in) :: x,p,phi0
    
    gdwi_epsilon_two = (8*p*(x*x)**p*(-1 + 2*p + (x*x)**p))/(phi0**2*x*x*(-1 + (x*x)**p)**2)

    
  end function gdwi_epsilon_two


!epsilon_three(x)
  function gdwi_epsilon_three(x,p,phi0)    
    implicit none
    real(kp) :: gdwi_epsilon_three
    real(kp), intent(in) :: x,p,phi0
    
    gdwi_epsilon_three = (8*p*(x*x)**p*(1 - 3*p + 2*p**2 + (2 + p)*(-1 + 2*p)*(x*x)**p &
         + (x*x)**(2*p)))/(phi0**2*x*x*(-1 + (x*x)**p)**2*(-1 + 2*p + (x*x)**p))
   
  end function gdwi_epsilon_three



  function gdwi_x_epsoneunity(p,phi0)
    implicit none
    real(kp), dimension(2) :: gdwi_x_epsoneunity
    real(kp), intent(in) :: p, phi0
    real(kp), parameter :: tolFind = tolkp
    
    real(kp) :: mini,maxi
    real(kp) :: xepstwonull

    type(transfert) :: gdwiData


    gdwiData%real1 = p
    gdwiData%real2 = phi0

    mini = 0._kp
    maxi = 1._kp - epsilon(1._kp)

    gdwi_x_epsoneunity(1) = zbrent(find_gdwi_x_epsoneunity,mini,maxi,tolFind,gdwiData)


    mini = 1._kp + epsilon(1._kp)
    maxi = huge(1._kp)

    gdwi_x_epsoneunity(2) = zbrent(find_gdwi_x_epsoneunity,mini,maxi,tolFind,gdwiData)

  end function gdwi_x_epsoneunity



  function find_gdwi_x_epsoneunity(x,gdwiData)
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: gdwiData
    real(kp) :: find_gdwi_x_epsoneunity
    real(kp) :: p, phi0

    p = gdwiData%real1
    phi0 = gdwiData%real2

    find_gdwi_x_epsoneunity = gdwi_epsilon_one(x,p,phi0) - 1._kp

  end function find_gdwi_x_epsoneunity
  


  function gdwi_x_endinf(p,phi0)
    implicit none
    real(kp), intent(in) :: p,phi0
    real(kp) :: gdwi_x_endinf
    real(kp), dimension(2) :: xepsone
    
    if (p.eq.1._kp) then
       gdwi_x_endinf = -sqrt(2._kp)/phi0+sqrt(1._kp+2._kp/phi0**2)
       return
    endif

    xepsone = gdwi_x_epsoneunity(p,phi0)
    
    gdwi_x_endinf = xepsone(1)
    

  end function gdwi_x_endinf




!this is integral[V(phi)/V'(phi) dphi]
  function gdwi_efold_primitive(x,p,phi0)
    implicit none
    real(kp), intent(in) :: x,p,phi0
    real(kp) :: gdwi_efold_primitive

    if (p.eq.1._kp) then
       gdwi_efold_primitive = 0.25_kp*phi0**2*(0.5_kp*x**2-log(x))
       return
    endif

    gdwi_efold_primitive =x*x*(-1+p+(x*x)**(-p))*phi0**2/(8*p*(p-1))
    
  end function gdwi_efold_primitive




  function gdwi_x_trajectory(bfold,xend,p,phi0)
    implicit none
    real(kp), intent(in) :: bfold,xend,p,phi0
    real(kp) :: gdwi_x_trajectory

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gdwiData

    real(kp), dimension(2) :: xepsone

    xepsone = gdwi_x_epsoneunity(p,phi0)

    mini = epsilon(0._kp)
    maxi = xepsone(1)

    gdwiData%real1 = p
    gdwiData%real2 = phi0
    gdwiData%real3 = -bfold + gdwi_efold_primitive(xend,p,phi0)

    gdwi_x_trajectory = zbrent(find_gdwi_x_trajectory,mini,maxi,tolFind,gdwiData)

  end function gdwi_x_trajectory

 


   function find_gdwi_x_trajectory(x,gdwiData)
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: gdwiData
    real(kp) :: find_gdwi_x_trajectory
    real(kp) :: p,phi0,NplusPrimEnd

    p = gdwiData%real1
    phi0 = gdwiData%real2
    NplusPrimEnd = gdwiData%real3

    find_gdwi_x_trajectory = gdwi_efold_primitive(x,p,phi0) - NplusPrimEnd

  end function find_gdwi_x_trajectory
  

end module gdwisr
