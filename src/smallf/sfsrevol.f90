!slow-roll functions for the small field potential
!
!V(phi) = M^4 [1 - (phi/mu)^p]
!
!x = phi/mu

module sfsrevol
  use infprec, only : kp, tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public  sf_norm_potential, sf_epsilon_one, sf_epsilon_two
  public  sf_x_endinf, sf_nufunc
 
contains
!returns V/M^4
  function sf_norm_potential(x,p)
    implicit none
    real(kp) :: sf_norm_potential
    real(kp), intent(in) :: x,p

    sf_norm_potential = 1._kp - x**p
  end function sf_norm_potential


!epsilon1(x)
  function sf_epsilon_one(x,p,mu)    
    implicit none
    real(kp) :: sf_epsilon_one
    real(kp), intent(in) :: x,p,mu
    
    sf_epsilon_one = 0.5_kp*(p/mu * (x**(p-1._kp))/(1._kp-x**p))**2
    
  end function sf_epsilon_one


!epsilon2(x)
  function sf_epsilon_two(x,p,mu)    
    implicit none
    real(kp) :: sf_epsilon_two
    real(kp), intent(in) :: x,p,mu
    
    sf_epsilon_two = 2._kp*(p/mu**2)*x**(p-2._kp) &
         * (p-1._kp + x**p)/(1._kp-x**p)**2
    
  end function sf_epsilon_two



!this is nu(x)=integral[V(phi)/V'(phi) dphi]
  function sf_nufunc(x,p,mu)
    implicit none
    real(kp), intent(in) :: x,p,mu
    real(kp) :: sf_nufunc
    
    if (p == 2._kp) then
       sf_nufunc = x**2 - 2._kp * log(x)
    else
       sf_nufunc = x**2 + 2._kp/(p-2._kp) * x**(2._kp-p)
    endif

    if (p.eq.0._kp) stop 'sf_nufunc: p=0 is singular'

    sf_nufunc = 0.5_kp*sf_nufunc*mu*mu/p

  end function sf_nufunc


  

!returns x at the end of inflation defined as epsilon1=1
  function sf_x_endinf(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: sf_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sfData

    mini = epsilon(1._kp)
    maxi = 1._kp + epsilon(1._kp)

    sfData%real1 = p
    sfData%real2 = mu

    sf_x_endinf = zbrent(find_sfendinf,mini,maxi,tolFind,sfData)
   
  end function sf_x_endinf
  
  function find_sfendinf(x,sfData)
    use infprec, only : transfert    
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: sfData
    real(kp) :: find_sfendinf
    real(kp) :: p,mu
    
    p=sfData%real1
    mu=sfData%real2
!this is epsilon1(x)=1
    find_sfendinf = x**(p-1._kp) + sqrt(2._kp)*mu/abs(p) * (x**p - 1._kp)
    
  end function find_sfendinf
 


!returns x at bfold=-efolds before the end of inflation
  function sf_x_trajectory(bfold,xend,p,mu)
    implicit none
    real(kp), intent(in) :: bfold, p, mu, xend
    real(kp) :: sf_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sfData

    mini = epsilon(1._kp)
    maxi = xend

    sfData%real1 = p
    sfData%real2 = mu
    sfData%real3 = -bfold + sf_nufunc(xend,p,mu)
    
    sf_x_trajectory = zbrent(find_sftraj,mini,maxi,tolFind,sfData)
    
   
  end function sf_x_trajectory

  function find_sftraj(x,sfData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: sfData
    real(kp) :: find_sftraj
    real(kp) :: p,mu,NplusNuend

    p=sfData%real1
    mu = sfData%real2
    NplusNuend = sfData%real3

    find_sftraj = sf_nufunc(x,p,mu) - NplusNuend
   
  end function find_sftraj

  
!returns x given epsilon2  
  function sf_x_epstwo(eps2,p,mu)   
    implicit none
    real(kp), intent(in) :: p,mu,eps2
    real(kp) :: sf_x_epstwo
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sfData

    mini = epsilon(1._kp)
    maxi = 1._kp + epsilon(1._kp)

    sfData%real1 = p
    sfData%real2 = mu
    sfData%real3 = eps2

   sf_x_epstwo = zbrent(find_sfepstwo,mini,maxi,tolFind,sfData)
   
 end function sf_x_epstwo

 function find_sfepstwo(x,sfData)    
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sfData
    real(kp) :: find_sfepstwo

    real(kp) :: p,mu,eps2

    p = sfData%real1
    mu = sfData%real2
    eps2 = sfData%real3

   find_sfepstwo  = sf_epsilon_two(x,p,mu) - eps2

 end function find_sfepstwo


end module sfsrevol
