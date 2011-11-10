!slow-roll functions for the mixed large field potential
!
!V(phi) = M^4 (phi/Mp)^p [1+alpha/q*(phi/Mp)^q]
!
!x = phi

module mlfisr
  use infprec, only : kp, tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : hypergeom_2F1
  implicit none

  private

  public  mlfi_norm_potential, mlfi_epsilon_one, mlfi_epsilon_two, mlfi_epsilon_three
  public  mlfi_x_endinf, mlfi_efold_primitive, mlfi_x_trajectory
 
contains
!returns V/M^4
  function mlfi_norm_potential(x,p,q,alpha)
    implicit none
    real(kp) :: mlfi_norm_potential
    real(kp), intent(in) :: x,p,q,alpha

    mlfi_norm_potential = x**p*(1._kp+alpha/q*x**q)

  end function mlfi_norm_potential


!epsilon_one(x)
  function mlfi_epsilon_one(x,p,q,alpha)    
    implicit none
    real(kp) :: mlfi_epsilon_one
    real(kp), intent(in) :: x,p,q,alpha
    
    mlfi_epsilon_one = 0.5_kp/(x**2) &
         *((p+alpha*(1._kp+p/q)*x**q)/(1._kp+alpha/q*x**q))**2
    
  end function mlfi_epsilon_one


!epsilon_two(x)
  function mlfi_epsilon_two(x,p,q,alpha)
    implicit none
    real(kp) :: mlfi_epsilon_two
    real(kp), intent(in) :: x,p,q,alpha
    
    mlfi_epsilon_two = 2_kp/(x**2) &
        *(p+alpha**2/(q**2)*(p+q)*x**(2._kp*q)+alpha/q*(2._kp*p+q-q**2)*x**q) &
         /(1._kp+alpha/q*x**q)**2  
    
  end function mlfi_epsilon_two


!epsilon_three(x)
  function mlfi_epsilon_three(x,p,q,alpha)
    implicit none
    real(kp) :: mlfi_epsilon_three
    real(kp), intent(in) :: x,p,q,alpha
    
    mlfi_epsilon_three = 1._kp/(x**2) &
         *1._kp/(p*q**2+alpha**2*(p+q)*x**(2*q)+alpha*q*(2._kp*p+q-q**2)*x**q) &
         *1._kp/(1._kp+alpha/q*x**q)**2 &
         *(2._kp*p**2*q**2+2._kp*alpha**4*(1+p/q)**2*x**(4._kp*q) &
         +alpha**2*(12._kp**2+6._kp*p*q*(2._kp-q)+(q-2._kp)*(q-1._kp)*q**2)*x**(2._kp*q) &
         +alpha**3*(p+q)*(8._kp*p/q+(1._kp-q)*(4._kp+q))*x**(3._kp*q) &
         +alpha*p*q*(8._kp*p+q*(4._kp+q**2-3._kp*q))*x**q)
    
  end function mlfi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function mlfi_x_endinf(p,q,alpha)
    implicit none
    real(kp), intent(in) :: p,q,alpha
    real(kp) :: mlfi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: mlfiData

    mini = epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    mlfiData%real1 = p
    mlfiData%real2 = q
    mlfiData%real3 = alpha

    mlfi_x_endinf = zbrent(find_mlfiendinf,mini,maxi,tolFind,mlfiData)
   
  end function mlfi_x_endinf


  function find_mlfiendinf(x,mlfiData)
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: mlfiData
    real(kp) :: find_mlfiendinf
    real(kp) :: p,q,alpha
    
    p = mlfiData%real1
    q = mlfiData%real2
    alpha = mlfiData%real3
    
    find_mlfiendinf = x*sqrt(2._kp)*(1._kp+alpha/q*x**q) &
         -(p+alpha*(1._kp+p/q)*x**q)
    
  end function find_mlfiendinf


!this is integral[V(phi)/V'(phi) dphi]
  function mlfi_efold_primitive(x,p,q,alpha)
    implicit none
    real(kp), intent(in) :: x,p,q,alpha
    real(kp) :: mlfi_efold_primitive
    
    mlfi_efold_primitive = 0.5_kp/(p+q)*x**2*(1._kp+q/p &
         *hypergeom_2F1(1._kp,2._kp/q,1._kp+2._kp/q,-alpha*(1._kp/p+1._kp/q)*x**q))
         
   
  end function mlfi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function mlfi_x_trajectory(bfold,xend,p,q,alpha)
    implicit none
    real(kp), intent(in) :: bfold, p, q, alpha, xend
    real(kp) :: mlfi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: mlfiData

  
    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)
  
    
    mlfiData%real1 = p
    mlfiData%real2 = q
    mlfiData%real3 = alpha
    mlfiData%real4 = -bfold + mlfi_efold_primitive(xend,p,q,alpha)
    
    mlfi_x_trajectory = zbrent(find_mlfitraj,mini,maxi,tolFind,mlfiData)
       
  end function mlfi_x_trajectory

  function find_mlfitraj(x,mlfiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: mlfiData
    real(kp) :: find_mlfitraj
    real(kp) :: p,q,alpha,NplusNuend

    p=mlfiData%real1
    q = mlfiData%real2
    alpha = mlfiData%real3
    NplusNuend = mlfiData%real4

    find_mlfitraj = mlfi_efold_primitive(x,p,q,alpha) - NplusNuend
   
  end function find_mlfitraj

  
end module mlfisr
