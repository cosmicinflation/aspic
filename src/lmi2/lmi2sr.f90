!slow-roll functions for the Logamediate inflation 2 potential ("2"
!means that inflation proceeds from the left to the right)
!
!V(phi) = M^4 x^[4(1-gamma)] exp[-beta * x^gamma]
!
!x = (phi-phi0)/Mp


module lmi2sr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : hypergeom_2F1
  use inftools, only : zbrent
  implicit none

  private

  public  lmi2_norm_potential, lmi2_epsilon_one, lmi2_epsilon_two, lmi2_epsilon_three
  public  lmi2_efold_primitive, lmi2_x_trajectory
  public  lmi2_norm_deriv_potential, lmi2_norm_deriv_second_potential
  public  lmi2_xin_min, lmi2_epsilon_one_max, lmi2_x_max_potential
 
contains

!returns V/M^4
  function lmi2_norm_potential(x,gamma,beta)
    implicit none
    real(kp) :: lmi2_norm_potential
    real(kp), intent(in) :: x,gamma,beta

    lmi2_norm_potential = x**(4.*(1.-gamma))*exp(-beta*x**gamma)

  end function lmi2_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function lmi2_norm_deriv_potential(x,gamma,beta)
    implicit none
    real(kp) :: lmi2_norm_deriv_potential
    real(kp), intent(in) :: x,gamma,beta
    real(kp) :: alpha

    alpha = 4._kp*(1._kp-gamma)

    lmi2_norm_deriv_potential = (alpha*x**(alpha-1._kp)- &
                                beta*gamma*x**(alpha+gamma-1._kp))*exp(-beta*x**gamma)
    
  end function lmi2_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function lmi2_norm_deriv_second_potential(x,gamma,beta)
    implicit none
    real(kp) :: lmi2_norm_deriv_second_potential
    real(kp), intent(in) :: x,gamma,beta
    real(kp) :: alpha

    alpha = 4._kp*(1._kp-gamma)

    lmi2_norm_deriv_second_potential = (x**(-2 + alpha)*((-1 + alpha)*alpha &
         - beta*gamma*(-1 + 2*alpha + gamma)*x**gamma &
         + beta**2*gamma**2*x**(2*gamma)))*exp(-beta*x**gamma)

  end function lmi2_norm_deriv_second_potential



!epsilon_one(x)
  function lmi2_epsilon_one(x,gamma,beta)    
    implicit none
    real(kp) :: lmi2_epsilon_one
    real(kp), intent(in) :: x,gamma,beta

    real(kp) ::alpha
    alpha=4.*(1.-gamma)

    lmi2_epsilon_one = (alpha-beta*gamma*x**gamma)**2/(2._kp*x**2)
    
  end function lmi2_epsilon_one


!epsilon_two(x)
  function lmi2_epsilon_two(x,gamma,beta)    
    implicit none
    real(kp) :: lmi2_epsilon_two
    real(kp), intent(in) :: x,gamma,beta

    real(kp) ::alpha
    alpha=4.*(1.-gamma)
    
    lmi2_epsilon_two = 2._kp*(alpha+beta*(-1._kp+gamma) &
                       *gamma*x**gamma)/(x**2)
    
  end function lmi2_epsilon_two


!epsilon_three(x)
  function lmi2_epsilon_three(x,gamma,beta)    
    implicit none
    real(kp) :: lmi2_epsilon_three
    real(kp), intent(in) :: x,gamma,beta

    real(kp) ::alpha
    alpha=4.*(1.-gamma)
    
    lmi2_epsilon_three = ((alpha-beta*gamma*x**gamma)*(2._kp*alpha- & 
                          beta*(-2._kp+gamma)*(-1._kp+gamma)*gamma* &
                          x**gamma))/(x**2*(alpha+beta*(-1._kp+gamma)* & 
                          gamma*x**gamma))
    
  end function lmi2_epsilon_three


  function lmi2_epsilon_one_max(gamma,beta)
    implicit none
    real(kp) :: lmi2_epsilon_one_max
    real(kp), intent(in) :: gamma,beta
    real(kp) :: alpha

    alpha = 4._kp*(1._kp - gamma)

    lmi2_epsilon_one_max= 0.5_kp*(alpha*gamma/(1._kp-gamma))**2 &
         *(beta*gamma*(1._kp-gamma)/alpha)**(2._kp/gamma)

  end function lmi2_epsilon_one_max


  function lmi2_x_max_potential(gamma,beta)
    implicit none
    real(kp) :: lmi2_x_max_potential
    real(kp), intent(in) :: gamma,beta
    real(kp) :: alpha

    alpha = 4._kp*(1._kp-gamma)

    lmi2_x_max_potential = (alpha/(beta*gamma))**(1._kp/gamma)

  end function lmi2_x_max_potential



!returns the minimum value for xin: if eps1<1 in the whole x>xVmax interval, returns xVmax, otherwise, returns the highest solution for eps1=1
  function lmi2_xin_min(gamma,beta)
    implicit none
    real(kp), intent(in) :: gamma,beta
    real(kp) :: lmi2_xin_min
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi2Data

    real(kp) ::alpha
    alpha=4._kp*(1._kp-gamma)

    if(beta.lt.sqrt(2._kp).or.gamma.lt.lmi2_gammamin(beta)) then
      lmi2_xin_min=(alpha/(beta*gamma))**(1._kp/gamma) &
                   *(1._kp+epsilon(1._kp))
    
    else

    mini = ((alpha/(beta*gamma*(1._kp-gamma)))**(1./gamma))*(1._kp+epsilon(1._kp))
    maxi = 100._kp*max(alpha,(beta*gamma)**(1._kp/(1._kp-gamma)),(alpha*beta*gamma)**(1._kp/(2._kp-gamma)))

    lmi2Data%real1 = gamma
    lmi2Data%real2 = beta
    
    lmi2_xin_min = zbrent(find_lmi2xinmin,mini,maxi,tolFind,lmi2Data)*(1._kp+epsilon(1._kp))

    endif


  end function lmi2_xin_min



  function find_lmi2xinmin(x,lmi2Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi2Data
    real(kp) :: find_lmi2xinmin
    real(kp) :: gamma,beta

    gamma = lmi2Data%real1
    beta = lmi2Data%real2
    
    find_lmi2xinmin = lmi2_epsilon_one(x,gamma,beta) - 1._kp
   
  end function find_lmi2xinmin




!this is integral[V(phi)/V'(phi) dphi]
  function lmi2_efold_primitive(x,gamma,beta)
    implicit none
    real(kp), intent(in) :: x,gamma,beta
    real(kp) :: lmi2_efold_primitive

    real(kp) ::alpha
    alpha=4._kp*(1._kp-gamma)

    if (alpha.eq.0._kp) stop 'lmi2_efold_primitive: gamma=1!  (PLI)'
    if (gamma.eq.0._kp) stop 'lmi2_efold_primitive: gamma=0!'

    lmi2_efold_primitive = x**2/(2._kp*alpha)*hypergeom_2F1(1._kp,2._kp/gamma, &
                           2._kp/gamma+1._kp,beta*gamma/alpha*x**gamma)


  end function lmi2_efold_primitive

!this is integral[V(phi)/V'(phi) dphi], approximated in the limit x/x0>>1
  function lmi2_efold_primitive_approximated(x,gamma,beta)
    implicit none
    real(kp), intent(in) :: x,gamma,beta
    real(kp) :: lmi2_efold_primitive_approximated

    real(kp) ::alpha
    alpha=4._kp*(1._kp-gamma)

    if (gamma.eq.0._kp) stop 'lmi2_efold_primitive: gamma=0!'

    lmi2_efold_primitive_approximated = x**(2._kp-gamma)/(beta*gamma*(gamma-2._kp))

  end function lmi2_efold_primitive_approximated


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function lmi2_x_trajectory(bfold,xend,gamma,beta)
    implicit none
    real(kp), intent(in) :: bfold, gamma, xend,beta
    real(kp) :: lmi2_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi, xiniMin
    type(transfert) :: lmi2Data

    real(kp) ::alpha
    alpha=4._kp*(1._kp-gamma)
    xiniMin = lmi2_xin_min(gamma,beta)

    if (xend.le.xiniMin) then
       write(*,*)'xiniMin= xend= ',xiniMin, xend
       stop 'lmi2_x_trajectory: xend < xiniMin'
    endif
    maxi = xend*(1._kp-epsilon(1._kp))
    mini = xiniMin

    lmi2Data%real1 = gamma
    lmi2Data%real2 = beta
    lmi2Data%real3 = -bfold + lmi2_efold_primitive(xend,gamma,beta)
    
    lmi2_x_trajectory = zbrent(find_lmi2traj,mini,maxi,tolFind,lmi2Data)
       
  end function lmi2_x_trajectory

  function find_lmi2traj(x,lmi2Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi2Data
    real(kp) :: find_lmi2traj
    real(kp) :: gamma,beta,NplusNuend

    gamma = lmi2Data%real1
    beta = lmi2Data%real2
    NplusNuend = lmi2Data%real3

    find_lmi2traj = lmi2_efold_primitive(x,gamma,beta) - NplusNuend

  end function find_lmi2traj


! Returns the minimum value for beta in order to end inflation with
! slow roll violation ( beta>betamin(gamma) <=> epsOneMax>1 )
  function lmi2_betamin(gamma)
    implicit none
    real(kp), intent(in) :: gamma
    real(kp) :: lmi2_betamin
    real(kp) ::alpha
    alpha=4._kp*(1._kp-gamma)

    lmi2_betamin=(sqrt(2._kp)*(1._kp-gamma)/(alpha*gamma))**gamma*alpha/(gamma*(1._kp-gamma))

  end function lmi2_betamin

! Returns the minimum value for gamma in order to end inflation with
! slow roll violation ( gamma>gammamin(beta) <=> epsOneMax>1 )
  function lmi2_gammamin(beta)
    implicit none
    real(kp), intent(in) :: beta
    real(kp) :: lmi2_gammamin
    real(kp), parameter :: tolFind=tolkp
    type(transfert) :: lmi2Data

    if (beta.lt.sqrt(2._kp)) stop 'lmi2_gammamin: beta<sqrt(2): inflation cannot end by slow-roll violation!'

    lmi2Data%real1 = beta

    lmi2_gammamin=zbrent(find_lmi2_gammamin,epsilon(1._kp),1._kp-epsilon(1._kp),tolFind,lmi2Data)

  end function lmi2_gammamin

  function find_lmi2_gammamin(x,lmi2Data)
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lmi2Data
    real(kp) :: find_lmi2_gammamin
    real(kp)  ::  beta

    beta = lmi2Data%real1
    
    find_lmi2_gammamin = lmi2_betamin(x)-beta

  end function find_lmi2_gammamin


end module lmi2sr
