!slow-roll functions for the Colemann Weinberg inflation potential
!
!V(phi) = M^4 * [1 + alpha x^4 ln(x) ]
!
!x = phi/Q
!alpha = 4e in order to set the potential minimum to 0.

module cwisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : ei,lambert
  use inftools, only : zbrent
  implicit none

  private

  public  cwi_norm_potential, cwi_epsilon_one, cwi_epsilon_two, cwi_epsilon_three
  public  cwi_x_endinf, cwi_efold_primitive, cwi_x_trajectory
  public  cwi_norm_deriv_potential, cwi_norm_deriv_second_potential
 
contains
!returns V/M^4
  function cwi_norm_potential(x,alpha,Q)
    implicit none
    real(kp) :: cwi_norm_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in), optional :: Q

    cwi_norm_potential = 1._kp+alpha*x**4._kp*log(x)

  end function cwi_norm_potential


![ cwi_xminus_positive_potential , cwi_xplus_positive_potential ] is
!the interval for x in which the potential is negative, if alpha is
!not set to its usual value

!Returns the lower bound of the x-interval within which the potential is negative
function cwi_xminus_positive_potential(alpha,Q)
    implicit none
    real(kp) :: cwi_xminus_positive_potential
    real(kp), intent(in) :: alpha
    real(kp), intent(in), optional :: Q

    if(alpha .lt. 4._kp*exp(1._kp) ) then
         cwi_xminus_positive_potential=1._kp/epsilon(1._kp)
    else
         cwi_xminus_positive_potential=exp(0.25_kp*lambert(-4._kp/(alpha),-1))
    endif

end function cwi_xminus_positive_potential

!Returns the upper bound of the x-interval within which the potential is negative
function cwi_xplus_positive_potential(alpha,Q)
    implicit none
    real(kp) :: cwi_xplus_positive_potential
    real(kp), intent(in) :: alpha
    real(kp), intent(in), optional :: Q

    if(alpha .lt. 4._kp*exp(1._kp) ) then
         cwi_xplus_positive_potential=0._kp
    else
         cwi_xplus_positive_potential=exp(0.25_kp*lambert(-4._kp/(alpha),0))
    endif

end function cwi_xplus_positive_potential


!returns the first derivative of the potential with respect to x=phi/Q, divided by M^4
  function cwi_norm_deriv_potential(x,alpha,Q)
    implicit none
    real(kp) :: cwi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in), optional :: Q

   cwi_norm_deriv_potential = alpha*x**3*(1._kp+4._kp*log(x))

  end function cwi_norm_deriv_potential



!returns the second derivative of the potential with respect to x=phi/Q, divided by M^4
  function cwi_norm_deriv_second_potential(x,alpha,Q)
    implicit none
    real(kp) :: cwi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in), optional :: Q

    cwi_norm_deriv_second_potential = alpha*x**2*(7._kp+12._kp*log(x))

  end function cwi_norm_deriv_second_potential



!epsilon_one(x)
  function cwi_epsilon_one(x,alpha,Q)    
    implicit none
    real(kp) :: cwi_epsilon_one
    real(kp), intent(in) :: x,alpha,Q

  
    cwi_epsilon_one = alpha**2/(2._kp*Q**2)*x**6* &
          ((1._kp+4._kp*log(x))/(1._kp+alpha*x**4*log(x)))**2

    
  end function cwi_epsilon_one


!epsilon_two(x)
  function cwi_epsilon_two(x,alpha,Q)    
    implicit none
    real(kp) :: cwi_epsilon_two
    real(kp), intent(in) :: x,alpha,Q
    
    cwi_epsilon_two = 2._kp*alpha/(Q**2)*x**2 &
         *(1._kp+alpha*x**4*log(x))**(-2) &
         *(-7._kp-12._kp*log(x)+alpha*x**4+alpha*x**4*log(x) &
         +4._kp*alpha*x**4*(log(x))**2)
    
  end function cwi_epsilon_two


!epsilon_three(x)
  function cwi_epsilon_three(x,alpha,Q)    
    implicit none
    real(kp) :: cwi_epsilon_three
    real(kp), intent(in) :: x,alpha,Q
    
    cwi_epsilon_three = 1._kp/(Q**2)*(-26._kp*alpha*x**2+21._kp*alpha**2*x**6 &
         -2._kp*alpha**3*x**10-128._kp*alpha*x**2*log(x) &
         +152._kp*alpha**2*x**6*log(x)-11*alpha**3*x**10*log(x) &
         -96._kp*alpha*x**2*(log(x))**2+368._kp*alpha**2*x**6*(log(x))**2 &
         -14._kp*alpha**3*x**10*(log(x))**2+384._kp*alpha**2*x**6*(log(x))**3 &
         -16._kp*alpha**3*x**10*(log(x))**3-32._kp*alpha**3*x**10*(log(x))**4) &
         *(1._kp+alpha*x**4*log(x))**(-2)*(7._kp-alpha*x**4+12._kp*log(x) &
         -alpha*x**4*log(x)-4._kp*alpha*x**4*(log(x))**2)**(-1)
    
  end function cwi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function cwi_x_endinf(alpha,Q)
    implicit none
    real(kp), intent(in) :: alpha,Q
    real(kp) :: cwi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: cwiData  

    maxi = exp(-1._kp/4._kp)*(1._kp-1000._kp*epsilon(1._kp))
    mini = epsilon(1._kp)*maxi

    cwiData%real1 = alpha
    cwiData%real2 = Q	
    
    cwi_x_endinf = zbrent(find_cwi_x_endinf,mini,maxi,tolFind,cwiData)

  end function cwi_x_endinf



  function find_cwi_x_endinf(x,cwiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: cwiData
    real(kp) :: find_cwi_x_endinf
    real(kp) :: alpha,Q

    alpha = cwiData%real1
    Q = cwiData%real2
    
    find_cwi_x_endinf = cwi_epsilon_one(x,alpha,Q) - 1._kp
   
  end function find_cwi_x_endinf




!this is integral[V(phi)/V'(phi) dphi]
  function cwi_efold_primitive(x,alpha,Q)
    implicit none
    real(kp), intent(in) :: x,alpha,Q
    real(kp) :: cwi_efold_primitive

    if (alpha.eq.0._kp) stop 'cwi_efold_primitive: alpha=0!'

    cwi_efold_primitive = Q**2*(sqrt(exp(1._kp))/(4._kp*alpha) &
                          *ei(-0.5_kp-2._kp*log(x)) &
                          -1._kp/(16._kp*sqrt(exp(1._kp))) &
                          *ei(0.5+2._kp*log(x))+0.125_kp*x**2)

  end function cwi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function cwi_x_trajectory(bfold,xend,alpha,Q)
    implicit none
    real(kp), intent(in) :: bfold, alpha,Q, xend
    real(kp) :: cwi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: cwiData

  
    mini = epsilon(1._kp)
    maxi = cwi_x_endinf(alpha,Q)
  

    cwiData%real1 = alpha
    cwiData%real2 = Q	
    cwiData%real3 = -bfold + cwi_efold_primitive(xend,alpha,Q)
    
    cwi_x_trajectory = zbrent(find_cwi_x_trajectory,mini,maxi,tolFind,cwiData)
       
  end function cwi_x_trajectory

  function find_cwi_x_trajectory(x,cwiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: cwiData
    real(kp) :: find_cwi_x_trajectory
    real(kp) :: alpha,Q,NplusNuend

    alpha = cwiData%real1
    Q = cwiData%real2
    NplusNuend = cwiData%real3

    find_cwi_x_trajectory = cwi_efold_primitive(x,alpha,Q) - NplusNuend
   
  end function find_cwi_x_trajectory



end module cwisr
