!slow-roll functions for the KKLT potential in the minimal kinetic
!term approximation
!
!V(phi) = M^4 / [1 + x^(-p)]
!
!x = phi/mu

module kkltsrevol
  use infprec, only : kp, tolkp,transfert
  use inftools, only : zbrent
  implicit none

!if true, find endinf with eps1=1, otherwise eps2=1
  logical, parameter :: endinfIsEpsOne=.true.

  private

  public kklt_norm_potential, kklt_epsilon_one, kklt_epsilon_two
  public kklt_x_endinf, kklt_nufunc, kklt_x_trajectory
  public kklt_x_epstwo
 
contains
!returns V/M^4
  function kklt_norm_potential(x,p)
    implicit none
    real(kp) :: kklt_norm_potential
    real(kp), intent(in) :: x,p

    kklt_norm_potential = 1._kp/(1._kp - x**(-p))

  end function kklt_norm_potential


!epsilon1(x)
  function kklt_epsilon_one(x,p,mu)    
    implicit none
    real(kp) :: kklt_epsilon_one
    real(kp), intent(in) :: x,p,mu
    
    kklt_epsilon_one = 0.5_kp*(p/mu /(x**(p+1._kp) + x))**2
    
  end function kklt_epsilon_one


!epsilon2(x)
  function kklt_epsilon_two(x,p,mu)    
    implicit none
    real(kp) :: kklt_epsilon_two
    real(kp), intent(in) :: x,p,mu
    
    kklt_epsilon_two = 2._kp*(p/mu**2)*x**(-p-2._kp) &
         * (p+1._kp + x**(-p)) &
         / (1._kp + x**(-p))**2
    
  end function kklt_epsilon_two



!this is nu(x)=integral[V(phi)/V'(phi) dphi]
  function kklt_nufunc(x,p,mu)
    implicit none
    real(kp), intent(in) :: x,p,mu
    real(kp) :: kklt_nufunc
    
    kklt_nufunc = 0.5_kp*x**2 + x**(p+2._kp)/(p+2._kp)
    kklt_nufunc = mu*mu * kklt_nufunc / p

  end function kklt_nufunc


  
!returns x at the end of inflation defined as epsilon1=1 or epsilon2=1
  function kklt_x_endinf(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: kklt_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kkltData

    mini = epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    kkltData%real1 = p
    kkltData%real2 = mu

    kklt_x_endinf = zbrent(find_kkltendinf,mini,maxi,tolFind,kkltData)
   
  end function kklt_x_endinf
  
  function find_kkltendinf(x,kkltData)
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: kkltData
    real(kp) :: find_kkltendinf
    real(kp) :: p,mu
    
    p=kkltData%real1
    mu=kkltData%real2

    if (endinfIsEpsOne) then
!eps1=1
       find_kkltendinf = x**(p+1._kp) + x - p/(mu*sqrt(2._kp))
    else
!eps2=1
       find_kkltendinf = (x**(p+1._kp) + x)**2 - (2._kp*p/mu/mu) &
            *( (p+1._kp)*x**p + 1._kp)
    endif
    
  end function find_kkltendinf
 


!returns x at bfold=-efolds before the end of inflation
  function kklt_x_trajectory(bfold,xend,p,mu)
    implicit none
    real(kp), intent(in) :: bfold, p, mu, xend
    real(kp) :: kklt_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kkltData

    mini = xend
    maxi = 1._kp/epsilon(1._kp)
    
    kkltData%real1 = p
    kkltData%real2 = mu
    kkltData%real3 = -bfold + kklt_nufunc(xend,p,mu)
    
    kklt_x_trajectory = zbrent(find_kklt_x_trajectory,mini,maxi,tolFind,kkltData)
       
  end function kklt_x_trajectory

  function find_kklt_x_trajectory(x,kkltData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: kkltData
    real(kp) :: find_kklt_x_trajectory
    real(kp) :: p,mu,NplusNuend

    p=kkltData%real1
    mu = kkltData%real2
    NplusNuend = kkltData%real3

    find_kklt_x_trajectory = kklt_nufunc(x,p,mu) - NplusNuend
   
  end function find_kklt_x_trajectory

  
!returns x given epsilon2  
  function kklt_x_epstwo(eps2,p,mu)   
    implicit none
    real(kp), intent(in) :: p,mu,eps2
    real(kp) :: kklt_x_epstwo
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kkltData

    mini = 1._kp + epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    kkltData%real1 = p
    kkltData%real2 = mu
    kkltData%real3 = eps2

    kklt_x_epstwo = zbrent(find_kkltepstwo,mini,maxi,tolFind,kkltData)
   
  end function kklt_x_epstwo

 function find_kkltepstwo(x,kkltData)    
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: kkltData
    real(kp) :: find_kkltepstwo

    real(kp) :: p,mu,eps2

    p = kkltData%real1
    mu = kkltData%real2
    eps2 = kkltData%real3

    find_kkltepstwo  = kklt_epsilon_two(x,p,mu) - eps2

  end function find_kkltepstwo

end module kkltsrevol
