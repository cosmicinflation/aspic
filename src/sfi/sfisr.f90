!slow-roll functions for the small field potential
!
!V(phi) = M^4 [1 - x^p]
!
!x = phi/mu
!mu = mu/Mp

module sfisr
  use infprec, only : kp, tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public sfi_norm_potential, sfi_norm_deriv_potential, sfi_norm_deriv_second_potential
  public sfi_epsilon_one, sfi_epsilon_two,sfi_epsilon_three
  public sfi_x_endinf, sfi_efold_primitive, sfi_x_trajectory

  public sfi_x_fromepstwo, find_sfi_x_fromepstwo, find_sfi_x_trajectory
  public find_sfi_x_endinf
 
contains
!returns V/M^4
  function sfi_norm_potential(x,p,mu)
    implicit none
    real(kp) :: sfi_norm_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in), optional :: mu

    sfi_norm_potential = 1._kp - x**p
  end function sfi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function sfi_norm_deriv_potential(x,p,mu)
    implicit none
    real(kp) :: sfi_norm_deriv_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in), optional :: mu

   sfi_norm_deriv_potential = -p*x**(p-1._kp)

  end function sfi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function sfi_norm_deriv_second_potential(x,p,mu)
    implicit none
    real(kp) :: sfi_norm_deriv_second_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in), optional :: mu

    sfi_norm_deriv_second_potential = -p*(p-1._kp)*x**(p-2._kp)

  end function sfi_norm_deriv_second_potential

!epsilon1(x)
  function sfi_epsilon_one(x,p,mu)    
    implicit none
    real(kp) :: sfi_epsilon_one
    real(kp), intent(in) :: x,p,mu
    
    sfi_epsilon_one = 0.5_kp*(p/mu * (x**(p-1._kp))/(1._kp-x**p))**2
    
  end function sfi_epsilon_one


!epsilon2(x)
  function sfi_epsilon_two(x,p,mu)    
    implicit none
    real(kp) :: sfi_epsilon_two
    real(kp), intent(in) :: x,p,mu
    
    sfi_epsilon_two = 2._kp*(p/mu**2)*x**(p-2._kp) &
         * (p-1._kp + x**p)/(1._kp-x**p)**2
    
  end function sfi_epsilon_two

!epsilon3(x)
  function sfi_epsilon_three(x,p,mu)    
    implicit none
    real(kp) :: sfi_epsilon_three
    real(kp), intent(in) :: x,p,mu
    
    sfi_epsilon_three = p/mu**2*x**(p-2)*(2*x**(2*p) &
    +(p-1)*(p+4)*x**p+(p-1)*(p-2))/((1-x**p)**2*(x**p+p-1))
    
  end function sfi_epsilon_three


!this is integral[V(phi)/V'(phi) dphi]
  function sfi_efold_primitive(x,p,mu)
    implicit none
    real(kp), intent(in) :: x,p,mu
    real(kp) :: sfi_efold_primitive
    
    if (p == 2._kp) then
       sfi_efold_primitive = x**2 - 2._kp * log(x)
    else
       sfi_efold_primitive = x**2 + 2._kp/(p-2._kp) * x**(2._kp-p)
    endif

    if (p.eq.0._kp) stop 'sf_nufunc: p=0 is singular'

    sfi_efold_primitive = 0.5_kp*sfi_efold_primitive*mu*mu/p

  end function sfi_efold_primitive


!returns x at the end of inflation defined as epsilon1=1
  function sfi_x_endinf(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: sfi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sfiData

    mini = epsilon(1._kp)
    maxi = 1._kp + epsilon(1._kp)

    sfiData%real1 = p
    sfiData%real2 = mu

    sfi_x_endinf = zbrent(find_sfi_x_endinf,mini,maxi,tolFind,sfiData)
   
  end function sfi_x_endinf
  
  function find_sfi_x_endinf(x,sfiData)
    use infprec, only : transfert    
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: sfiData
    real(kp) :: find_sfi_x_endinf
    real(kp) :: p,mu
    
    p=sfiData%real1
    mu=sfiData%real2
!this is epsilon1(x)=1
    find_sfi_x_endinf = x**(p-1._kp) + sqrt(2._kp)*mu/abs(p) * (x**p - 1._kp)
    
  end function find_sfi_x_endinf
 


!returns x at bfold=-efolds before the end of inflation
  function sfi_x_trajectory(bfold,xend,p,mu)
    implicit none
    real(kp), intent(in) :: bfold, p, mu, xend
    real(kp) :: sfi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sfiData

    mini = epsilon(1._kp)
    maxi = xend

    sfiData%real1 = p
    sfiData%real2 = mu
    sfiData%real3 = -bfold + sfi_efold_primitive(xend,p,mu)
    
    sfi_x_trajectory = zbrent(find_sfi_x_trajectory,mini,maxi,tolFind,sfiData)
    
  end function sfi_x_trajectory

  function find_sfi_x_trajectory(x,sfiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: sfiData
    real(kp) :: find_sfi_x_trajectory
    real(kp) :: p,mu,NplusPrimEnd

    p=sfiData%real1
    mu = sfiData%real2
    NplusPrimEnd = sfiData%real3

    find_sfi_x_trajectory = sfi_efold_primitive(x,p,mu) - NplusPrimEnd
   
  end function find_sfi_x_trajectory

  
!returns x given epsilon2  
  function sfi_x_fromepstwo(eps2,p,mu)   
    implicit none
    real(kp), intent(in) :: p,mu,eps2
    real(kp) :: sfi_x_fromepstwo
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sfiData

    mini = epsilon(1._kp)
    maxi = 1._kp + epsilon(1._kp)

    sfiData%real1 = p
    sfiData%real2 = mu
    sfiData%real3 = eps2

   sfi_x_fromepstwo = zbrent(find_sfi_x_fromepstwo,mini,maxi,tolFind,sfiData)
   
 end function sfi_x_fromepstwo

 function find_sfi_x_fromepstwo(x,sfiData)    
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sfiData
    real(kp) :: find_sfi_x_fromepstwo

    real(kp) :: p,mu,eps2

    p = sfiData%real1
    mu = sfiData%real2
    eps2 = sfiData%real3

   find_sfi_x_fromepstwo  = sfi_epsilon_two(x,p,mu) - eps2

 end function find_sfi_x_fromepstwo


end module sfisr
