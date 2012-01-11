!slow-roll functions for the loop inflation potential
!
!V(phi) = M^4 [1 + alpha ln(x) ]
!
!x = phi/Mp

module lisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  use inftools, only : zbrent
  implicit none

  private

  public  li_norm_potential, li_epsilon_one, li_epsilon_two, li_epsilon_three
  public  li_x_endinf, li_efold_primitive, li_x_trajectory
  public  li_norm_deriv_potential, li_norm_deriv_second_potential
 
contains
!returns V/M^4
  function li_norm_potential(x,alpha)
    implicit none
    real(kp) :: li_norm_potential
    real(kp), intent(in) :: x,alpha

    li_norm_potential = 1._kp+alpha*log(x)

  end function li_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function li_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: li_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

   li_norm_deriv_potential = alpha/x

  end function li_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function li_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: li_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    li_norm_deriv_second_potential = -alpha/(x**2)

  end function li_norm_deriv_second_potential



!epsilon_one(x)
  function li_epsilon_one(x,alpha)
    implicit none
    real(kp) :: li_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    li_epsilon_one = alpha**2/(2._kp*x**2) &
         /(1._kp+alpha*log(x))**2
    
  end function li_epsilon_one


!epsilon_two(x)
  function li_epsilon_two(x,alpha)
    implicit none
    real(kp) :: li_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    li_epsilon_two = 2._kp*alpha/(x**2)*(1._kp+alpha+alpha*log(x)) &
         /(1._kp+alpha*log(x))**2
    
  end function li_epsilon_two


!epsilon_three(x)
  function li_epsilon_three(x,alpha)
    implicit none
    real(kp) :: li_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    li_epsilon_three = 2._kp*alpha/(x**2) &
         /(1._kp+alpha+alpha*log(x)) &
         /(1._kp+alpha*log(x))**2 &
         *(1._kp+3._kp*alpha/2._kp+alpha**2+ &
         (2._kp*alpha+3._kp*alpha**2/2._kp)*log(x) &
         +alpha**2*(log(x))**2)
    
  end function li_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function li_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: li_x_endinf
    

    li_x_endinf = 1._kp/sqrt(2._kp) &
         /lambert(exp(1._kp/alpha)/(sqrt(2._kp)),0)

   
  end function li_x_endinf


  

!this is integral[V(phi)/V'(phi) dphi]
  function li_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: li_efold_primitive

    if (alpha.eq.0._kp) stop 'li_efold_primitive: alpha=0!'

    li_efold_primitive = (-1._kp/4._kp+1._kp/(2._kp*alpha))*x**2 &
         +1._kp/2._kp*x**2*log(x)

  end function li_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function li_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: li_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: liData

  
    mini = xEnd
    maxi = xEnd*1000._kp
  


    liData%real1 = alpha
    liData%real2 = -bfold + li_efold_primitive(xend,alpha)
    
    li_x_trajectory = zbrent(find_litraj,mini,maxi,tolFind,liData)
       
  end function li_x_trajectory

  function find_litraj(x,liData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: liData
    real(kp) :: find_litraj
    real(kp) :: alpha,mu,NplusNuend

    alpha = liData%real1
    NplusNuend = liData%real2

    find_litraj = li_efold_primitive(x,alpha) - NplusNuend
   
  end function find_litraj


  
end module lisr
