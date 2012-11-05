!slow-roll functions for the generalized mixed large field inflation potential
!
!V(phi) = M^4 x^p [ 1 + alpha x^q ]
!
!x = phi/Mp

module gmlfisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : hypergeom_2F1

  implicit none

  private

  public gmlfi_norm_potential, gmlfi_epsilon_one, gmlfi_epsilon_two
  public gmlfi_epsilon_three
  public gmlfi_x_endinf, gmlfi_efold_primitive, gmlfi_x_trajectory
  public gmlfi_norm_deriv_potential, gmlfi_norm_deriv_second_potential

 
contains
!returns V/M**4
  function gmlfi_norm_potential(x,p,q,alpha)
    implicit none
    real(kp) :: gmlfi_norm_potential
    real(kp), intent(in) :: x,alpha,p,q

    gmlfi_norm_potential = x**p*(1._kp+alpha*x**q)

  end function gmlfi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function gmlfi_norm_deriv_potential(x,p,q,alpha)
    implicit none
    real(kp) :: gmlfi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,p,q

   gmlfi_norm_deriv_potential = x**(p-1._kp)*(p+alpha*(p+q)*x**q)

  end function gmlfi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function gmlfi_norm_deriv_second_potential(x,p,q,alpha)
    implicit none
    real(kp) :: gmlfi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,p,q

    gmlfi_norm_deriv_second_potential = x**(p-2._kp)*(p*(p-1._kp) &
         +alpha*(p+q)*(p+q-1._kp)*x**q)

  end function gmlfi_norm_deriv_second_potential



!epsilon_one(x)
  function gmlfi_epsilon_one(x,p,q,alpha)    
    implicit none
    real(kp) :: gmlfi_epsilon_one
    real(kp), intent(in) :: x,alpha,p,q
    
    gmlfi_epsilon_one = (p+alpha*(p+q)*x**q)**2 &
         /(2._kp*x**2*(1._kp+alpha*x**q)**2)
    
  end function gmlfi_epsilon_one


!epsilon_two(x)
  function gmlfi_epsilon_two(x,p,q,alpha)    
    implicit none
    real(kp) :: gmlfi_epsilon_two
    real(kp), intent(in) :: x,alpha,p,q
    
    gmlfi_epsilon_two = (2._kp*(p*(1._kp+alpha*x**q)**2 &
         +alpha*q*x**q*(1._kp-q+alpha*x**q))) & 
         /(x**2*(1._kp+alpha*x**q)**2)
    
  end function gmlfi_epsilon_two


!epsilon_three(x)
  function gmlfi_epsilon_three(x,p,q,alpha)    
    implicit none
    real(kp) :: gmlfi_epsilon_three
    real(kp), intent(in) :: x,alpha,p,q
    
    gmlfi_epsilon_three = ((p+alpha*(p+q)*x**q)*(2._kp*p*(1._kp+alpha*x**q)**3+ &
         alpha*q*x**q*(2._kp-3._kp*q+q**2-alpha*(-1._kp+q)*(4._kp+q)*x**q &
         +2._kp*alpha**2*x**(2._kp*q))))/(x**2*(1._kp+alpha*x**q)**2* &
         (p*(1._kp+alpha*x**q)**2+alpha*q*x**q*(1._kp-q+alpha*x**q)))

    
  end function gmlfi_epsilon_three


!returns the value for x=phi/Mp defined as epsilon1=1, where inflation ends
  function gmlfi_x_endinf(p,q,alpha)
    implicit none
    real(kp), intent(in) :: alpha,p,q
    real(kp) :: gmlfi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gmlfiData

    mini = p/10._kp**(6._kp) ! since eps1 scales as p²/x² when x->0
    maxi = (p+q)*10._kp**(6._kp) ! since eps1 scales as (p+q)²/x² when x->0 number
  

    gmlfiData%real1 = alpha
    gmlfiData%real2 = p
    gmlfiData%real3 = q
    
    gmlfi_x_endinf = zbrent(find_gmlfi_x_endinf,mini,maxi,tolFind,gmlfiData)
   
   
  end function gmlfi_x_endinf

  function find_gmlfi_x_endinf(x,gmlfiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: gmlfiData
    real(kp) :: find_gmlfi_x_endinf
    real(kp) :: alpha,p,q

    alpha = gmlfiData%real1
    p = gmlfiData%real2
    q = gmlfiData%real3

    find_gmlfi_x_endinf = gmlfi_epsilon_one(x,p,q,alpha)-1._kp
   
  end function find_gmlfi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function gmlfi_efold_primitive(x,p,q,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha,p,q
    real(kp) :: gmlfi_efold_primitive

    gmlfi_efold_primitive = 0.5_kp/(p+q)*x**2*(1+q/p* &
         hypergeom_2F1(1._kp,2._kp/q,1._kp+2._kp/q,-alpha*q*(1._kp/p+1._kp/q)*x**q))

  end function gmlfi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function gmlfi_x_trajectory(bfold,xend,p,q,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha,p,q, xend
    real(kp) :: gmlfi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gmlfiData

  
    mini = xend
    maxi = xend*10._kp**(6._kp)

    gmlfiData%real1 = alpha
    gmlfiData%real2 = p
    gmlfiData%real3 = q
    gmlfiData%real4 = -bfold + gmlfi_efold_primitive(xend,p,q,alpha)
    
    gmlfi_x_trajectory = zbrent(find_gmlfitraj,mini,maxi,tolFind,gmlfiData)
       
  end function gmlfi_x_trajectory

  function find_gmlfitraj(x,gmlfiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: gmlfiData
    real(kp) :: find_gmlfitraj
    real(kp) :: alpha,p,q,NplusNuend

    alpha= gmlfiData%real1
    p= gmlfiData%real2
    q= gmlfiData%real3
    NplusNuend = gmlfiData%real4

    find_gmlfitraj = gmlfi_efold_primitive(x,p,q,alpha) - NplusNuend
   
  end function find_gmlfitraj

  
end module gmlfisr
