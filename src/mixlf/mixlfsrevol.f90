!slow-roll functions for the large field mixed potential
!
!V(phi) = M^4 [phi^p + alpha phi^q]
!
!x = phi

module mixlfsrevol
  use infprec, only : kp, tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : hypergeom_2F1
  implicit none

  private

  public mixlf_norm_potential, mixlf_epsilon_one, mixlf_epsilon_two
  public mixlf_x_endinf, mixlf_nufunc, mixlf_x_trajectory
 
contains
!returns V/M^4
  function mixlf_norm_potential(x,p,q,alpha)
    implicit none
    real(kp) :: mixlf_norm_potential
    real(kp), intent(in) :: x,p,q,alpha

    mixlf_norm_potential = x**p + alpha * x**q

  end function mixlf_norm_potential


!epsilon1(x)
  function mixlf_epsilon_one(x,p,q,alpha)    
    implicit none
    real(kp) :: mixlf_epsilon_one
    real(kp), intent(in) :: x,p,q,alpha
    
    mixlf_epsilon_one = 0.5_kp/x**2 &
         *((p+alpha*q*x**(q-p))/(1+alpha*x**(q-p)))**2
    
  end function mixlf_epsilon_one


!epsilon2(x)
  function mixlf_epsilon_two(x,p,q,alpha)
    implicit none
    real(kp) :: mixlf_epsilon_two
    real(kp), intent(in) :: x,p,q,alpha
    
    mixlf_epsilon_two = 2_kp/x**2 &
         * (p-alpha*(p*p+q*(q-1)-p*(2*q+1))*x**(q-p) &
         + alpha*alpha*q*(x*x)**(q-p)) &
         / (1+alpha*x**(q-p))**2
    
  end function mixlf_epsilon_two



!this is nu(x)=integral[V(phi)/V'(phi) dphi]
  function mixlf_nufunc(x,p,q,alpha)
    implicit none
    real(kp), intent(in) :: x,p,q,alpha
    real(kp) :: mixlf_nufunc
    
    mixlf_nufunc = 0.5_kp*x**2/(p*q)*(p + (q-p) &
         *hypergeom_2F1(1._kp,2._kp/(q-p),1._kp + 2._kp/(q-p),-alpha*q*x**(q-p)/p))
   
  end function mixlf_nufunc


  
!returns x at the end of inflation defined as epsilon1=1 or epsilon2=1
  function mixlf_x_endinf(p,q,alpha)
    implicit none
    real(kp), intent(in) :: p,q,alpha
    real(kp) :: mixlf_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: mixlfData

    mini = epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    mixlfData%real1 = p
    mixlfData%real2 = q
    mixlfData%real3 = alpha

    mixlf_x_endinf = zbrent(find_mixlfendinf,mini,maxi,tolFind,mixlfData)
   
  end function mixlf_x_endinf
  
  function find_mixlfendinf(x,mixlfData)
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: mixlfData
    real(kp) :: find_mixlfendinf
    real(kp) :: p,q,alpha
    
    p = mixlfData%real1
    q = mixlfData%real2
    alpha = mixlfData%real3
    
    find_mixlfendinf = alpha*sqrt(2._kp)*x**(q-p+1._kp) - alpha*q*x**(q-p) &
         + sqrt(2._kp)*x - p
    
  end function find_mixlfendinf
 


!returns x at bfold=-efolds before the end of inflation
  function mixlf_x_trajectory(bfold,xend,p,q,alpha)
    implicit none
    real(kp), intent(in) :: bfold, p, q, alpha, xend
    real(kp) :: mixlf_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: mixlfData

  
    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)
  
    
    mixlfData%real1 = p
    mixlfData%real2 = q
    mixlfData%real3 = alpha
    mixlfData%real4 = -bfold + mixlf_nufunc(xend,p,q,alpha)
    
    mixlf_x_trajectory = zbrent(find_mixlftraj,mini,maxi,tolFind,mixlfData)
       
  end function mixlf_x_trajectory

  function find_mixlftraj(x,mixlfData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: mixlfData
    real(kp) :: find_mixlftraj
    real(kp) :: p,q,alpha,NplusNuend

    p=mixlfData%real1
    q = mixlfData%real2
    alpha = mixlfData%real3
    NplusNuend = mixlfData%real4

    find_mixlftraj = mixlf_nufunc(x,p,q,alpha) - NplusNuend
   
  end function find_mixlftraj

  
end module mixlfsrevol
