!slow-roll functions for the Kähler moduli inflation I potential
!
!V(phi) = M^4 [1 - alpha x exp(-x) ]
!
!x = phi/Mp

module kmiisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert,ei
  use inftools, only : zbrent
  implicit none

  private

  public  kmii_norm_potential, kmii_epsilon_one, kmii_epsilon_two, kmii_epsilon_three
  public  kmii_x_endinf, kmii_efold_primitive, kmii_x_trajectory
  public  kmii_norm_deriv_potential, kmii_norm_deriv_second_potential
 
contains
!returns V/M^4
  function kmii_norm_potential(x,alpha)
    implicit none
    real(kp) :: kmii_norm_potential
    real(kp), intent(in) :: x,alpha

    kmii_norm_potential = 1._kp-alpha*x*exp(-x)

  end function kmii_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function kmii_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) ::  kmii_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

   kmii_norm_deriv_potential = -alpha*exp(-x)+alpha*x*exp(-x)

  end function kmii_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function kmii_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: kmii_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    kmii_norm_deriv_second_potential = 2._kp*alpha*exp(-x)-alpha*x*exp(-x)

  end function kmii_norm_deriv_second_potential



!epsilon_one(x)
  function kmii_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: kmii_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    kmii_epsilon_one = alpha**2/2._kp * exp(-2._kp*x) &
         *(1._kp-x)**2/(1._kp-alpha*x*exp(-x))**2
    
  end function kmii_epsilon_one


!epsilon_two(x)
  function kmii_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: kmii_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    kmii_epsilon_two = 2._kp*alpha*exp(-x) &
         *(alpha*exp(-x)+x-2._kp) &
         /(1._kp-alpha*x*exp(-x))**2
    
  end function kmii_epsilon_two


!epsilon_three(x)
  function kmii_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: kmii_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    kmii_epsilon_three = alpha*exp(-x)*(x-1._kp) &
         /(1._kp-alpha*x*exp(-x))**2 &
         /(alpha*exp(-x)+x-2._kp) &
         *(x-3._kp+alpha*exp(-x)*(x**2-3._kp*x+6._kp) &
         -2._kp*alpha**2*exp(-2._kp*x))
    
  end function kmii_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function kmii_x_endinf(alpha) !assuming that inflation proceeds from the right to the left,
				!from initial high values of the field compared with the Planck mass
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: kmii_x_endinf
    real(kp) :: alpha1,alpha2 !delimitation if the different regimes
    
    alpha1=sqrt(2._kp)/(sqrt(2._kp)-1._kp) &
         *exp((2._kp-sqrt(2._kp))/(1._kp-sqrt(2._kp)))

    alpha2=sqrt(2._kp)/(sqrt(2._kp)+1._kp) &
         *exp((2._kp+sqrt(2._kp))/(1._kp+sqrt(2._kp)))
    
    if (alpha.le.0) then
        print*, 'alpha has to be positive'

    else if (alpha.le.alpha1) then
         print*, 'alpha<alpha1=',alpha1,': Never Ending Inflation Regime. Need to add instability.'

    else if (alpha.le.alpha2) then
         print*, 'alpha1<alpha<alpha2',': Never Ending Inflation Regime &
                    (if inflation proceeds from the right to the left), &
                    or non slow-roll Inflation Regime (if inflation proceeds from the left to the right)'

    else if (alpha.le.exp(1._kp)) then
    	kmii_x_endinf = 1._kp/(1._kp+sqrt(2._kp)) &
             -lambert(-sqrt(2._kp)/(1._kp+sqrt(2._kp)) &
   	     *exp(1._kp/(1._kp+sqrt(2._kp)))/alpha,-1)

    else
        print*, 'alpha>e',':Non positive potential'

    end if
   
  end function kmii_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function kmii_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: kmii_efold_primitive

    if (alpha.eq.0._kp) stop 'kmii_efold_primitive: alpha=0!'

    kmii_efold_primitive = -x+exp(1._kp)/alpha*ei(x-1._kp) &
         -log(x-1._kp)

  end function kmii_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function kmii_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: kmii_x_trajectory
    real(kp), parameter :: tokmiind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kmiiData

  
    mini = xEnd
    maxi = 1._kp/epsilon(1._kp)
  
    kmiiData%real1 = alpha
    kmiiData%real2 = -bfold + kmii_efold_primitive(xend,alpha)
    
    kmii_x_trajectory = zbrent(find_kmiitraj,mini,maxi,tokmiind,kmiiData)
       
  end function kmii_x_trajectory

  function find_kmiitraj(x,kmiiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: kmiiData
    real(kp) :: find_kmiitraj
    real(kp) :: alpha,NplusNuend

    alpha = kmiiData%real1
    NplusNuend = kmiiData%real2

    find_kmiitraj = kmii_efold_primitive(x,alpha) - NplusNuend
   
  end function find_kmiitraj


  
end module kmiisr
