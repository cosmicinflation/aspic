!slow-roll functions for the generalized mixed inflation potential
!
!V(phi) = M^4 x^p [ 1 + alpha x^q ]
!
!x = phi/Mp

module gmisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : hypergeom_2F1

  implicit none

  private

  public  gmi_norm_potential, gmi_epsilon_one, gmi_epsilon_two, gmi_epsilon_three
  public  gmi_x_endinf, gmi_efold_primitive, gmi_x_trajectory
  public  gmi_norm_deriv_potential, gmi_norm_deriv_second_potential

 
contains
!returns V/M**4
  function gmi_norm_potential(x,alpha,p,q)
    implicit none
    real(kp) :: gmi_norm_potential
    real(kp), intent(in) :: x,alpha,p,q

    gmi_norm_potential = x**p*(1._kp+alpha*x**q)

  end function gmi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function gmi_norm_deriv_potential(x,alpha,p,q)
    implicit none
    real(kp) :: gmi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,p,q

   gmi_norm_deriv_potential = x**(p-1._kp)*(p+alpha*(p+q)*x**q)

  end function gmi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function gmi_norm_deriv_second_potential(x,alpha,p,q)
    implicit none
    real(kp) :: gmi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,p,q

    gmi_norm_deriv_second_potential = x**(p-2._kp)*(p*(p-1._kp)+alpha*(p+q)*(p+q-1._kp)*x**q)

  end function gmi_norm_deriv_second_potential



!epsilon_one(x)
  function gmi_epsilon_one(x,alpha,p,q)    
    implicit none
    real(kp) :: gmi_epsilon_one
    real(kp), intent(in) :: x,alpha,p,q
    
    gmi_epsilon_one = (p+alpha*(p+q)*x**q)**2/(2._kp*x**2*(1._kp+alpha*x**q)**2)
    
  end function gmi_epsilon_one


!epsilon_two(x)
  function gmi_epsilon_two(x,alpha,p,q)    
    implicit none
    real(kp) :: gmi_epsilon_two
    real(kp), intent(in) :: x,alpha,p,q
    
    gmi_epsilon_two = (2._kp*(p*(1._kp+alpha*x**q)**2+alpha*q*x**q*(1._kp-q+alpha*x**q))) & 
                      /(x**2*(1._kp+alpha*x**q)**2)
    
  end function gmi_epsilon_two


!epsilon_three(x)
  function gmi_epsilon_three(x,alpha,p,q)    
    implicit none
    real(kp) :: gmi_epsilon_three
    real(kp), intent(in) :: x,alpha,p,q
    
    gmi_epsilon_three = ((p+alpha*(p+q)*x**q)*(2._kp*p*(1._kp+alpha*x**q)**3+ &
                        alpha*q*x**q*(2._kp-3._kp*q+q**2-alpha*(-1._kp+q)*(4._kp+q)*x**q &
                        +2._kp*alpha**2*x**(2._kp*q))))/(x**2*(1._kp+alpha*x**q)**2* &
                        (p*(1._kp+alpha*x**q)**2+alpha*q*x**q*(1._kp-q+alpha*x**q)))

    
  end function gmi_epsilon_three


!returns the value for x=phi/Mp defined as epsilon1=1, where inflation ends
  function gmi_x_endinf(alpha,p,q)
    implicit none
    real(kp), intent(in) :: alpha,p,q
    real(kp) :: gmi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gmiData

    mini = p/10._kp**(6._kp) ! since eps1 scales as p²/x² when x->0
    maxi = (p+q)*10._kp**(6._kp) ! since eps1 scales as (p+q)²/x² when x->0 number
  

    gmiData%real1 = alpha
    gmiData%real2 = p
    gmiData%real3 = q
    
    gmi_x_endinf = zbrent(find_gmi_x_endinf,mini,maxi,tolFind,gmiData)
   
   
  end function gmi_x_endinf

  function find_gmi_x_endinf(x,gmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: gmiData
    real(kp) :: find_gmi_x_endinf
    real(kp) :: alpha,p,q

    alpha = gmiData%real1
    p = gmiData%real2
    q = gmiData%real3

    find_gmi_x_endinf = gmi_epsilon_one(x,alpha,p,q)-1._kp
   
  end function find_gmi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function gmi_efold_primitive(x,alpha,p,q)
    implicit none
    real(kp), intent(in) :: x,alpha,p,q
    real(kp) :: gmi_efold_primitive

    gmi_efold_primitive = 0.5_kp/(p+q)*x**2*(1+q/p* &
                          hypergeom_2F1(1._kp,2._kp/q,1._kp+2._kp/q,-alpha*q*(1._kp/p+1._kp/q)*x**q))

  end function gmi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function gmi_x_trajectory(bfold,xend,alpha,p,q)
    implicit none
    real(kp), intent(in) :: bfold, alpha,p,q, xend
    real(kp) :: gmi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gmiData

  
    mini = xend
    maxi = xend*10._kp**(6._kp)

    gmiData%real1 = alpha
    gmiData%real2 = p
    gmiData%real3 = q
    gmiData%real4 = -bfold + gmi_efold_primitive(xend,alpha,p,q)
    
    gmi_x_trajectory = zbrent(find_gmitraj,mini,maxi,tolFind,gmiData)
       
  end function gmi_x_trajectory

  function find_gmitraj(x,gmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: gmiData
    real(kp) :: find_gmitraj
    real(kp) :: alpha,p,q,NplusNuend

    alpha= gmiData%real1
    p= gmiData%real2
    q= gmiData%real3
    NplusNuend = gmiData%real4

    find_gmitraj = gmi_efold_primitive(x,alpha,p,q) - NplusNuend
   
  end function find_gmitraj

  
end module gmisr
