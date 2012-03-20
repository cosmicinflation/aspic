!slow-roll functions for the Logamediate inflation potential, 2nd regime
!
!V(phi) = M^4 * (phi/Mp)^alpha * exp[-beta * (phi/Mp)^gamma_lmi ]
!
!x = phi/Mp


module lmi2sr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : hypergeom_2F1
  use inftools, only : zbrent
  implicit none

  private

  public  lmi2_norm_potential, lmi2_epsilon_one, lmi2_epsilon_two, lmi2_epsilon_three
  public  lmi2_x_endinf, lmi2_efold_primitive, lmi2_x_trajectory
  public  lmi2_norm_deriv_potential, lmi2_norm_deriv_second_potential
  public  lmi2_beta,lmi2_alpha
 
contains

!Build the model parameters from beta consistently with the model ansatz a(t)=exp[A(ln(t)^lambda)]
  function lmi2_alpha(gamma_lmi,M)
    implicit none
    real(kp) :: lmi2_alpha
    real(kp), intent(in) :: gamma_lmi,M

    lmi2_alpha=4._kp*(1._kp-gamma_lmi/2._kp)

  end function lmi2_alpha

  function lmi2_beta(gamma_lmi,M)
    implicit none
    real(kp) :: lmi2_beta
    real(kp), intent(in) :: gamma_lmi,M

    lmi2_beta=sqrt(2._kp*sqrt(3._kp))/(M*gamma_lmi)

  end function lmi2_beta


!returns V/M^4
  function lmi2_norm_potential(x,gamma_lmi,M)
    implicit none
    real(kp) :: lmi2_norm_potential
    real(kp), intent(in) :: x,gamma_lmi,M

    real(kp) ::alpha,beta
    alpha=lmi2_alpha(gamma_lmi,M)
    beta=lmi2_beta(gamma_lmi,M)

    lmi2_norm_potential = x**alpha*exp(-beta*x**gamma_lmi)

  end function lmi2_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function lmi2_norm_deriv_potential(x,gamma_lmi,M)
    implicit none
    real(kp) :: lmi2_norm_deriv_potential
    real(kp), intent(in) :: x,gamma_lmi,M

    real(kp) ::alpha,beta
    alpha=lmi2_alpha(gamma_lmi,M)
    beta=lmi2_beta(gamma_lmi,M)

    lmi2_norm_deriv_potential = exp(-beta*x**gamma_lmi)*x**(-1._kp+alpha)* &
                               (alpha-beta*gamma_lmi*x**gamma_lmi)

  end function lmi2_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function lmi2_norm_deriv_second_potential(x,gamma_lmi,M)
    implicit none
    real(kp) :: lmi2_norm_deriv_second_potential
    real(kp), intent(in) :: x,gamma_lmi,M

    real(kp) ::alpha,beta
    alpha=lmi2_alpha(gamma_lmi,M)
    beta=lmi2_beta(gamma_lmi,M)

    lmi2_norm_deriv_second_potential = exp(-beta*x**gamma_lmi)*x**(-2._kp &
                                       +alpha)*(alpha**2-alpha*(1._kp+2._kp* &
                                       beta*gamma_lmi*x**gamma_lmi)+beta*gamma_lmi* &
                                       x**gamma_lmi*(1._kp+gamma_lmi*(-1._kp+beta*x**gamma_lmi)))

  end function lmi2_norm_deriv_second_potential



!epsilon_one(x)
  function lmi2_epsilon_one(x,gamma_lmi,M)    
    implicit none
    real(kp) :: lmi2_epsilon_one
    real(kp), intent(in) :: x,gamma_lmi,M

    real(kp) ::alpha,beta
    alpha=lmi2_alpha(gamma_lmi,M)
    beta=lmi2_beta(gamma_lmi,M)

    lmi2_epsilon_one = (alpha-beta*gamma_lmi*x**gamma_lmi)**2/(2._kp*x**2)
    
  end function lmi2_epsilon_one


!epsilon_two(x)
  function lmi2_epsilon_two(x,gamma_lmi,M)    
    implicit none
    real(kp) :: lmi2_epsilon_two
    real(kp), intent(in) :: x,gamma_lmi,M

    real(kp) ::alpha,beta
    alpha=lmi2_alpha(gamma_lmi,M)
    beta=lmi2_beta(gamma_lmi,M)
    
    lmi2_epsilon_two = 2._kp*(alpha+beta*(-1._kp+gamma_lmi) &
                       *gamma_lmi*x**gamma_lmi)/(x**2)
    
  end function lmi2_epsilon_two


!epsilon_three(x)
  function lmi2_epsilon_three(x,gamma_lmi,M)    
    implicit none
    real(kp) :: lmi2_epsilon_three
    real(kp), intent(in) :: x,gamma_lmi,M

    real(kp) ::alpha,beta
    alpha=lmi2_alpha(gamma_lmi,M)
    beta=lmi2_beta(gamma_lmi,M)
    
    lmi2_epsilon_three = ((alpha-beta*gamma_lmi*x**gamma_lmi)*(2._kp*alpha- & 
                          beta*(-2._kp+gamma_lmi)*(-1._kp+gamma_lmi)*gamma_lmi* &
                          x**gamma_lmi))/(x**2*(alpha+beta*(-1._kp+gamma_lmi)* & 
                          gamma_lmi*x**gamma_lmi))
    
  end function lmi2_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function lmi2_x_endinf(gamma_lmi,M)
    implicit none
    real(kp), intent(in) :: gamma_lmi,M
    real(kp) :: lmi2_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi2Data

    real(kp) ::alpha,beta
    alpha=lmi2_alpha(gamma_lmi,M)
    beta=lmi2_beta(gamma_lmi,M)

    maxi = (alpha/(beta*gamma_lmi*(1._kp-gamma_lmi)))**(1._kp/gamma_lmi)*(1._kp-epsilon(1._kp))
    mini = (alpha/(beta*gamma_lmi))**(1._kp/gamma_lmi)*(1._kp+epsilon(1._kp))

    lmi2Data%real1 = gamma_lmi
    lmi2Data%real2 = M
    
    lmi2_x_endinf = zbrent(find_lmi2endinf,mini,maxi,tolFind,lmi2Data)


  end function lmi2_x_endinf



  function find_lmi2endinf(x,lmi2Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi2Data
    real(kp) :: find_lmi2endinf
    real(kp) :: gamma_lmi,M

    gamma_lmi = lmi2Data%real1
    M = lmi2Data%real2
    
    find_lmi2endinf = lmi2_epsilon_one(x,gamma_lmi,M) - 1._kp
   
  end function find_lmi2endinf




!this is integral[V(phi)/V'(phi) dphi]
  function lmi2_efold_primitive(x,gamma_lmi,M)
    implicit none
    real(kp), intent(in) :: x,gamma_lmi,M
    real(kp) :: lmi2_efold_primitive

    real(kp) ::alpha,beta
    alpha=lmi2_alpha(gamma_lmi,M)
    beta=lmi2_beta(gamma_lmi,M)

    if (alpha.eq.0._kp) stop 'lmi2_efold_primitive: alpha=0!'
    if (gamma_lmi.eq.0._kp) stop 'lmi2_efold_primitive: gamma=0!'

    lmi2_efold_primitive = x**2/(2._kp*alpha)*hypergeom_2F1(1._kp,2._kp/gamma_lmi, &
                           2._kp/(gamma_lmi+1._kp),beta*gamma_lmi/alpha*x**gamma_lmi)

  end function lmi2_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function lmi2_x_trajectory(bfold,xend,gamma_lmi,M)
    implicit none
    real(kp), intent(in) :: bfold, gamma_lmi, xend,M
    real(kp) :: lmi2_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi2Data

  
    real(kp) ::alpha,beta
    alpha=lmi2_alpha(gamma_lmi,M)
    beta=lmi2_beta(gamma_lmi,M)

    maxi = lmi2_x_endinf(gamma_lmi,M)*(1._kp-epsilon(1._kp))
    mini = (alpha/(beta*gamma_lmi))**(1._kp/gamma_lmi)*(1._kp+epsilon(1._kp))
  

    lmi2Data%real1 = gamma_lmi
    lmi2Data%real2 = M
    lmi2Data%real3 = -bfold + lmi2_efold_primitive(xend,gamma_lmi,M)
    
    lmi2_x_trajectory = zbrent(find_lmi2traj,mini,maxi,tolFind,lmi2Data)
       
  end function lmi2_x_trajectory

  function find_lmi2traj(x,lmi2Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi2Data
    real(kp) :: find_lmi2traj
    real(kp) :: gamma_lmi,M,NplusNuend

    gamma_lmi = lmi2Data%real1
    M = lmi2Data%real2
    NplusNuend = lmi2Data%real2

    find_lmi2traj = lmi2_efold_primitive(x,gamma_lmi,M) - NplusNuend
   
  end function find_lmi2traj



end module lmi2sr
