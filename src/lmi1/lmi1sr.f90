!slow-roll functions for the Logamediate inflation 1 potential ("1" means that inflation proceeds from the right to the left)
!
!V(phi) = M^4 * x^[4(1-gamma_lmi)] * exp[-beta * x^gamma_lmi ]
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
 
contains

!returns V/M^4
  function lmi1_norm_potential(x,gamma_lmi,beta)
    implicit none
    real(kp) :: lmi1_norm_potential
    real(kp), intent(in) :: x,gamma_lmi,beta

    lmi1_norm_potential = x**(4.*(1.-gamma_lmi))*exp(-beta*x**gamma_lmi)

  end function lmi1_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function lmi1_norm_deriv_potential(x,gamma_lmi,beta)
    implicit none
    real(kp) :: lmi1_norm_deriv_potential
    real(kp), intent(in) :: x,gamma_lmi,beta

    lmi1_norm_deriv_potential = (4.*(1.-gamma_lmi)*x**(3.-4.*gamma_lmi)- &
                                beta*gamma_lmi*x**(3.*(gamma_lmi-1.)))*exp(-beta*x**gamma_lmi)

  end function lmi1_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function lmi1_norm_deriv_second_potential(x,gamma_lmi,beta)
    implicit none
    real(kp) :: lmi1_norm_deriv_second_potential
    real(kp), intent(in) :: x,gamma_lmi,beta

    lmi1_norm_deriv_second_potential = (4.*(1.-gamma_lmi)*x**(2.-4.*gamma_lmi)+ &
                                       beta**2*gamma_lmi**2*x**(2.-2.*gamma_lmi)- &
                                       beta*gamma_lmi*(5.-4.*gamma_lmi)*x**(2.-3.*gamma_lmi)) &
                                       *exp(-beta*x**gamma_lmi)

  end function lmi1_norm_deriv_second_potential



!epsilon_one(x)
  function lmi1_epsilon_one(x,gamma_lmi,beta)    
    implicit none
    real(kp) :: lmi1_epsilon_one
    real(kp), intent(in) :: x,gamma_lmi,beta

    real(kp) ::alpha
    alpha=4.*(1.-gamma_lmi)

    lmi1_epsilon_one = (alpha-beta*gamma_lmi*x**gamma_lmi)**2/(2._kp*x**2)
    
  end function lmi1_epsilon_one


!epsilon_two(x)
  function lmi1_epsilon_two(x,gamma_lmi,beta)    
    implicit none
    real(kp) :: lmi1_epsilon_two
    real(kp), intent(in) :: x,gamma_lmi,beta

    real(kp) ::alpha
    alpha=4.*(1.-gamma_lmi)
    
    lmi1_epsilon_two = 2._kp*(alpha+beta*(-1._kp+gamma_lmi) &
                       *gamma_lmi*x**gamma_lmi)/(x**2)
    
  end function lmi1_epsilon_two


!epsilon_three(x)
  function lmi1_epsilon_three(x,gamma_lmi,beta)    
    implicit none
    real(kp) :: lmi1_epsilon_three
    real(kp), intent(in) :: x,gamma_lmi,beta

    real(kp) ::alpha
    alpha=4.*(1.-gamma_lmi)
    
    lmi1_epsilon_three = ((alpha-beta*gamma_lmi*x**gamma_lmi)*(2._kp*alpha- & 
                          beta*(-2._kp+gamma_lmi)*(-1._kp+gamma_lmi)*gamma_lmi* &
                          x**gamma_lmi))/(x**2*(alpha+beta*(-1._kp+gamma_lmi)* & 
                          gamma_lmi*x**gamma_lmi))
    
  end function lmi1_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function lmi1_x_endinf(gamma_lmi,beta)
    implicit none
    real(kp), intent(in) :: gamma_lmi,beta
    real(kp) :: lmi1_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi1Data

    real(kp) ::alpha
    alpha=4.*(1.-gamma_lmi)

    maxi = (alpha/(beta*gamma_lmi))**(1._kp/gamma_lmi)*(1._kp-epsilon(1._kp))
    mini = epsilon(1._kp)

    lmi1Data%real1 = gamma_lmi
    lmi1Data%real2 = beta
    
    lmi1_x_endinf = zbrent(find_lmi1endinf,mini,maxi,tolFind,lmi1Data)


  end function lmi1_x_endinf



  function find_lmi1endinf(x,lmi1Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi1Data
    real(kp) :: find_lmi1endinf
    real(kp) :: gamma_lmi,beta

    gamma_lmi = lmi1Data%real1
    beta = lmi1Data%real2
    
    find_lmi1endinf = lmi1_epsilon_one(x,gamma_lmi,beta) - 1._kp
   
  end function find_lmi1endinf




!this is integral[V(phi)/V'(phi) dphi]
  function lmi1_efold_primitive(x,gamma_lmi,beta)
    implicit none
    real(kp), intent(in) :: x,gamma_lmi,beta
    real(kp) :: lmi1_efold_primitive

    real(kp) ::alpha
    alpha=4.*(1.-gamma_lmi)

    if (alpha.eq.0._kp) stop 'lmi1_efold_primitive: gamma=1!  (PLI)'
    if (gamma_lmi.eq.0._kp) stop 'lmi1_efold_primitive: gamma=0!'

    lmi1_efold_primitive = x**2/(2._kp*alpha)*hypergeom_2F1(1._kp,2._kp/gamma_lmi, &
                           2._kp/gamma_lmi+1._kp,beta*gamma_lmi/alpha*x**gamma_lmi)

  end function lmi1_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function lmi1_x_trajectory(bfold,xend,gamma_lmi,beta)
    implicit none
    real(kp), intent(in) :: bfold, gamma_lmi, xend,beta
    real(kp) :: lmi1_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi1Data

    real(kp) ::alpha
    alpha=4.*(1.-gamma_lmi)

    maxi = (alpha/(beta*gamma_lmi))**(1._kp/gamma_lmi)*(1._kp-epsilon(1._kp))
    mini = lmi1_x_endinf(gamma_lmi,beta)*(1._kp+epsilon(1._kp))
  

    lmi1Data%real1 = gamma_lmi
    lmi1Data%real2 = beta
    lmi1Data%real3 = -bfold + lmi1_efold_primitive(xend,gamma_lmi,beta)
    
    lmi1_x_trajectory = zbrent(find_lmi1traj,mini,maxi,tolFind,lmi1Data)
       
  end function lmi1_x_trajectory

  function find_lmi1traj(x,lmi1Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi1Data
    real(kp) :: find_lmi1traj
    real(kp) :: gamma_lmi,beta,NplusNuend

    gamma_lmi = lmi1Data%real1
    beta = lmi1Data%real2
    NplusNuend = lmi1Data%real2

    find_lmi1traj = lmi1_efold_primitive(x,gamma_lmi,beta) - NplusNuend
   
  end function find_lmi1traj



end module lmi1sr
