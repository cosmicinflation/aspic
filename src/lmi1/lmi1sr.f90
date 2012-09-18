!slow-roll functions for the Logamediate inflation 1 potential ("1"
!means that inflation proceeds from the right to the left)
!
!V(phi) = M^4 x^[4(1-gamma)] exp[-beta x^gamma]
!
!x = (phi-phi0)/Mp


module lmi1sr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : hypergeom_2F1
  use inftools, only : zbrent
  implicit none

  private

  public  lmi1_norm_potential, lmi1_epsilon_one, lmi1_epsilon_two, lmi1_epsilon_three
  public  lmi1_x_endinf, lmi1_efold_primitive, lmi1_x_trajectory
  public  lmi1_norm_deriv_potential, lmi1_norm_deriv_second_potential
  public  lmi1_epsilon_one_max, lmi1_x_max_potential

contains

!returns V/M^4
  function lmi1_norm_potential(x,gamma,beta)
    implicit none
    real(kp) :: lmi1_norm_potential
    real(kp), intent(in) :: x,gamma,beta
    real(kp) :: alpha

    alpha = 4._kp*(1._kp-gamma)

    lmi1_norm_potential = x**alpha*exp(-beta*x**gamma)

  end function lmi1_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function lmi1_norm_deriv_potential(x,gamma,beta)
    implicit none
    real(kp) :: lmi1_norm_deriv_potential
    real(kp), intent(in) :: x,gamma,beta
    real(kp) :: alpha

    alpha = 4._kp*(1._kp-gamma)

    lmi1_norm_deriv_potential = (alpha*x**(alpha-1._kp)- &
                                beta*gamma*x**(alpha+gamma-1._kp))*exp(-beta*x**gamma)

  end function lmi1_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function lmi1_norm_deriv_second_potential(x,gamma,beta)
    implicit none
    real(kp) :: lmi1_norm_deriv_second_potential
    real(kp), intent(in) :: x,gamma,beta
    real(kp) :: alpha

    alpha = 4._kp*(1._kp-gamma)
   
    lmi1_norm_deriv_second_potential = (x**(-2 + alpha)*((-1 + alpha)*alpha &
         - beta*gamma*(-1 + 2*alpha + gamma)*x**gamma &
         + beta**2*gamma**2*x**(2*gamma)))*exp(-beta*x**gamma)
    

  end function lmi1_norm_deriv_second_potential



!epsilon_one(x)
  function lmi1_epsilon_one(x,gamma,beta)    
    implicit none
    real(kp) :: lmi1_epsilon_one
    real(kp), intent(in) :: x,gamma,beta

    real(kp) ::alpha
    alpha=4._kp*(1._kp-gamma)

    lmi1_epsilon_one = (alpha-beta*gamma*x**gamma)**2/(2._kp*x**2)
    
  end function lmi1_epsilon_one


!epsilon_two(x)
  function lmi1_epsilon_two(x,gamma,beta)    
    implicit none
    real(kp) :: lmi1_epsilon_two
    real(kp), intent(in) :: x,gamma,beta

    real(kp) ::alpha
    alpha=4._kp*(1._kp-gamma)
    
    lmi1_epsilon_two = 2._kp*(alpha+beta*(-1._kp+gamma) &
                       *gamma*x**gamma)/(x**2)
    
  end function lmi1_epsilon_two


!epsilon_three(x)
  function lmi1_epsilon_three(x,gamma,beta)    
    implicit none
    real(kp) :: lmi1_epsilon_three
    real(kp), intent(in) :: x,gamma,beta

    real(kp) ::alpha
    alpha=4._kp*(1._kp-gamma)
    
    lmi1_epsilon_three = ((alpha-beta*gamma*x**gamma)*(2._kp*alpha- & 
                          beta*(-2._kp+gamma)*(-1._kp+gamma)*gamma* &
                          x**gamma))/(x**2*(alpha+beta*(-1._kp+gamma)* & 
                          gamma*x**gamma))
    
  end function lmi1_epsilon_three


  function lmi1_epsilon_one_max(gamma,beta)
    implicit none
    real(kp) :: lmi1_epsilon_one_max
    real(kp), intent(in) :: gamma,beta
    real(kp) :: alpha

    alpha = 4._kp*(1._kp - gamma)

    lmi1_epsilon_one_max= 0.5_kp*(alpha*gamma/(1._kp-gamma))**2 &
         *(beta*gamma*(1._kp-gamma)/alpha)**(2._kp/gamma)

  end function lmi1_epsilon_one_max


  function lmi1_x_max_potential(gamma,beta)
    implicit none
    real(kp) :: lmi1_x_max_potential
    real(kp), intent(in) :: gamma,beta
    real(kp) :: alpha

    alpha = .4_kp*(1._kp-gamma)

    lmi1_x_max_potential = alpha/(beta*gamma)**(1._kp/gamma)

  end function lmi1_x_max_potential



!returns x at the end of inflation defined as epsilon1=1
  function lmi1_x_endinf(gamma,beta)
    implicit none
    real(kp), intent(in) :: gamma,beta
    real(kp) :: lmi1_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi1Data

    real(kp) ::alpha
    alpha=4._kp*(1._kp-gamma)

    maxi = (alpha/(beta*gamma))**(1._kp/gamma)*(1._kp-epsilon(1._kp))
    mini = epsilon(1._kp)

    lmi1Data%real1 = gamma
    lmi1Data%real2 = beta
    
    lmi1_x_endinf = zbrent(find_lmi1endinf,mini,maxi,tolFind,lmi1Data)


  end function lmi1_x_endinf



  function find_lmi1endinf(x,lmi1Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi1Data
    real(kp) :: find_lmi1endinf
    real(kp) :: gamma,beta

    gamma = lmi1Data%real1
    beta = lmi1Data%real2
    
    find_lmi1endinf = lmi1_epsilon_one(x,gamma,beta) - 1._kp
   
  end function find_lmi1endinf




!this is integral[V(phi)/V'(phi) dphi]
  function lmi1_efold_primitive(x,gamma,beta)
    implicit none
    real(kp), intent(in) :: x,gamma,beta
    real(kp) :: lmi1_efold_primitive

    real(kp) ::alpha
    alpha=4._kp*(1._kp-gamma)

    if (alpha.eq.0._kp) stop 'lmi1_efold_primitive: gamma=1!  (PLI)'
    if (gamma.eq.0._kp) stop 'lmi1_efold_primitive: gamma=0!'

    lmi1_efold_primitive = x**2/(2._kp*alpha)*hypergeom_2F1(1._kp,2._kp/gamma, &
                           2._kp/gamma+1._kp,beta*gamma/alpha*x**gamma)

  end function lmi1_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function lmi1_x_trajectory(bfold,xend,gamma,beta)
    implicit none
    real(kp), intent(in) :: bfold, gamma, xend,beta
    real(kp) :: lmi1_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi1Data

    real(kp) ::alpha
    alpha=4._kp*(1._kp-gamma)

    maxi = (alpha/(beta*gamma))**(1._kp/gamma)*(1._kp-epsilon(1._kp))
    mini = lmi1_x_endinf(gamma,beta)*(1._kp+epsilon(1._kp))
  

    lmi1Data%real1 = gamma
    lmi1Data%real2 = beta
    lmi1Data%real3 = -bfold + lmi1_efold_primitive(xend,gamma,beta)
    
    lmi1_x_trajectory = zbrent(find_lmi1traj,mini,maxi,tolFind,lmi1Data)
       
  end function lmi1_x_trajectory

  function find_lmi1traj(x,lmi1Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi1Data
    real(kp) :: find_lmi1traj
    real(kp) :: gamma,beta,NplusNuend

    gamma = lmi1Data%real1
    beta = lmi1Data%real2
    NplusNuend = lmi1Data%real3

    find_lmi1traj = lmi1_efold_primitive(x,gamma,beta) - NplusNuend
   
  end function find_lmi1traj



end module lmi1sr
