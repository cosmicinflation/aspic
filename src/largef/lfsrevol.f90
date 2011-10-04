!slow-roll functions for the large field potential
!
!V(phi) = M^4 phi^p
!
!x = phi

module lfsrevol
  use infprec, only : kp
  implicit none

  private

  public  lf_norm_potential, lf_epsilon_one, lf_epsilon_two
  public  lf_x_endinf, lf_nufunc, lf_x_trajectory
 
contains
!returns V/M^4
  function lf_norm_potential(x,p)
    implicit none
    real(kp) :: lf_norm_potential
    real(kp), intent(in) :: x,p

    lf_norm_potential = x**p

  end function lf_norm_potential


!epsilon1(x)
  function lf_epsilon_one(x,p)    
    implicit none
    real(kp) :: lf_epsilon_one
    real(kp), intent(in) :: x,p
    
    lf_epsilon_one = 0.5_kp*(p/x)**2
    
  end function lf_epsilon_one


!epsilon2(x)
  function lf_epsilon_two(x,p)    
    implicit none
    real(kp) :: lf_epsilon_two
    real(kp), intent(in) :: x,p
    
    lf_epsilon_two = 2._kp*p/x**2
    
  end function lf_epsilon_two



!this is nu(x)=integral[V(phi)/V'(phi) dphi]
  function lf_nufunc(x,p)
    implicit none
    real(kp), intent(in) :: x,p
    real(kp) :: lf_nufunc

    if (p.eq.0._kp) stop 'lf_nufunc: p=0!'

    lf_nufunc = 0.5_kp*x**2/p

  end function lf_nufunc


  

!returns x at the end of inflation defined as epsilon1=1
  function lf_x_endinf(p)
    implicit none
    real(kp), intent(in) :: p
    real(kp) :: lf_x_endinf
   
    lf_x_endinf = p/sqrt(2._kp)
   
  end function lf_x_endinf
  
  


!returns x at bfold=-efolds before the end of inflation
  function lf_x_trajectory(bfold,xend,p)
    implicit none
    real(kp), intent(in) :: bfold, p, xend
    real(kp) :: lf_x_trajectory
           
    lf_x_trajectory = sqrt(-2._kp*p*bfold + xend**2)
    
   
  end function lf_x_trajectory

  
  
end module lfsrevol
