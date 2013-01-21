!slow-roll functions for the radiatively corrected quartic inflation potential
!
!V(phi) = M^4 x^4*[1-alpha ln(x) ]
!
!x = phi/Mp

module rcqisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : ei
  use inftools, only : zbrent
  implicit none

  private

  public  rcqi_norm_potential, rcqi_epsilon_one, rcqi_epsilon_two, rcqi_epsilon_three
  public  rcqi_x_endinf, rcqi_efold_primitive, rcqi_x_trajectory
  public  rcqi_norm_deriv_potential, rcqi_norm_deriv_second_potential
 
contains
!returns V/M^4
  function rcqi_norm_potential(x,alpha)
    implicit none
    real(kp) :: rcqi_norm_potential
    real(kp), intent(in) :: x,alpha

    rcqi_norm_potential = x**4*(1._kp-alpha*log(x))

  end function rcqi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function rcqi_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: rcqi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

   rcqi_norm_deriv_potential = 4._kp*x**3*(1._kp-alpha*log(x)) &
         -alpha*x**3

  end function rcqi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function rcqi_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: rcqi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    rcqi_norm_deriv_second_potential = 12._kp*x**2*(1._kp-alpha*log(x)) &
         -7._kp*alpha*x**2

  end function rcqi_norm_deriv_second_potential



!epsilon_one(x)
  function rcqi_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: rcqi_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    rcqi_epsilon_one = 8._kp/(x**2) &
         *(1._kp-alpha/4._kp-alpha*log(x))**2 &
         /(1._kp-alpha*log(x))**2
    
  end function rcqi_epsilon_one


!epsilon_two(x)
  function rcqi_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: rcqi_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    rcqi_epsilon_two = 2._kp/(x**2) &
         /(1._kp-alpha*log(x))**2 &
         *(4._kp-alpha+alpha**2-8._kp*alpha*log(x) &
         +alpha**2*log(x)+4._kp*alpha**2*(log(x))**2)
    
  end function rcqi_epsilon_two


!epsilon_three(x)
  function rcqi_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: rcqi_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    rcqi_epsilon_three = 1._kp/(x**2) &
         *(4._kp-alpha-4._kp*alpha*log(x)) &
         /(1._kp-alpha*log(x))**2 &
         *(8._kp-2._kp*alpha+3._kp*alpha**2-2._kp*alpha**3-24._kp*alpha*log(x) &
         +4._kp*alpha**2*log(x)+24._kp*alpha**2*(log(x))**2-3._kp*alpha**3*log(x)-2._kp*alpha**3*(log(x))**2 &
         -8._kp*alpha**3*(log(x))**3) &
         /(4._kp-alpha+alpha**2-8._kp*alpha*log(x)+alpha**2*log(x)+4._kp*alpha**2*(log(x))**2)
    
  end function rcqi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function rcqi_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: rcqi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rcqiData

    mini = epsilon(1._kp)
    maxi = min(1._kp/epsilon(1._kp),exp(-0.25_kp+1._kp/alpha))

    rcqiData%real1 = alpha

    rcqi_x_endinf = zbrent(find_rcqiendinf,mini,maxi,tolFind,rcqiData)
   
  end function rcqi_x_endinf


  function find_rcqiendinf(x,rcqiData)
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: rcqiData
    real(kp) :: find_rcqiendinf
    real(kp) :: alpha
    
    alpha = rcqiData%real1
    
    find_rcqiendinf = 2._kp*sqrt(2._kp)/x &
         *(1._kp-alpha/4._kp-alpha*log(x)) &
         +alpha*log(x)-1._kp
    
  end function find_rcqiendinf


!this is integral[V(phi)/V'(phi) dphi]
  function rcqi_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: rcqi_efold_primitive

    if (alpha.eq.0._kp) stop 'rcqi_efold_primitive: alpha=0!'

    rcqi_efold_primitive = 1._kp/16._kp*(2._kp*x**2-exp(-0.5_kp+2._kp/alpha) &
         *ei(0.5_kp-2._kp/alpha+2._kp*log(x)))

  end function rcqi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function rcqi_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: rcqi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rcqiData

  
    mini = xEnd
    maxi = min(1._kp/epsilon(1._kp),exp(-0.25_kp+1._kp/alpha))
  


    rcqiData%real1 = alpha
    rcqiData%real2 = -bfold + rcqi_efold_primitive(xend,alpha)
    
    rcqi_x_trajectory = zbrent(find_rcqi_x_trajectory,mini,maxi,tolFind,rcqiData)
       
  end function rcqi_x_trajectory

  function find_rcqi_x_trajectory(x,rcqiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: rcqiData
    real(kp) :: find_rcqi_x_trajectory
    real(kp) :: alpha,NplusNuend

    alpha = rcqiData%real1
    NplusNuend = rcqiData%real2

    find_rcqi_x_trajectory = rcqi_efold_primitive(x,alpha) - NplusNuend
   
  end function find_rcqi_x_trajectory


  
end module rcqisr
