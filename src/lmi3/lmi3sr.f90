!slow-roll functions for the Logamediate inflation potential, third regime of inflation
!
!V(phi) = M^4 x^alpha * exp[ -beta x^gamma_lmi ]
!
!x = phi/Mp


module lmi3sr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : hypergeom_2F1
  use inftools, only : zbrent
  implicit none

  private

  public  lmi3_norm_potential, lmi3_epsilon_one, lmi3_epsilon_two, lmi3_epsilon_three
  public  lmi3_efold_primitive, lmi3_x_trajectory
  public  lmi3_norm_deriv_potential, lmi3_norm_deriv_second_potential
  public  lmi3_beta,lmi3_alpha
 
contains

!Build the model parameters from beta consistently with the model ansatz a(t)=exp[A(ln(t)^lambda)]
  function lmi3_alpha(gamma_lmi,M)
    implicit none
    real(kp) :: lmi3_alpha
    real(kp), intent(in) :: gamma_lmi,M

    lmi3_alpha=4._kp*(1._kp-gamma_lmi/2._kp)

  end function lmi3_alpha

  function lmi3_beta(gamma_lmi,M)
    implicit none
    real(kp) :: lmi3_beta
    real(kp), intent(in) :: gamma_lmi,M

    lmi3_beta=sqrt(2._kp*sqrt(3._kp))/(M*gamma_lmi)

  end function lmi3_beta


!returns V/M^4
  function lmi3_norm_potential(x,gamma_lmi,M)
    implicit none
    real(kp) :: lmi3_norm_potential
    real(kp), intent(in) :: x,gamma_lmi,M

    real(kp) ::alpha,beta
    alpha=lmi3_alpha(gamma_lmi,M)
    beta=lmi3_beta(gamma_lmi,M)

    lmi3_norm_potential = x**alpha*exp(-beta*x**gamma_lmi)

  end function lmi3_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function lmi3_norm_deriv_potential(x,gamma_lmi,M)
    implicit none
    real(kp) :: lmi3_norm_deriv_potential
    real(kp), intent(in) :: x,gamma_lmi,M

    real(kp) ::alpha,beta
    alpha=lmi3_alpha(gamma_lmi,M)
    beta=lmi3_beta(gamma_lmi,M)

    lmi3_norm_deriv_potential = exp(-beta*x**gamma_lmi)*x**(-1._kp+alpha)* &
                               (alpha-beta*gamma_lmi*x**gamma_lmi)

  end function lmi3_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function lmi3_norm_deriv_second_potential(x,gamma_lmi,M)
    implicit none
    real(kp) :: lmi3_norm_deriv_second_potential
    real(kp), intent(in) :: x,gamma_lmi,M

    real(kp) ::alpha,beta
    alpha=lmi3_alpha(gamma_lmi,M)
    beta=lmi3_beta(gamma_lmi,M)

    lmi3_norm_deriv_second_potential = exp(-beta*x**gamma_lmi)*x**(-2._kp &
                                       +alpha)*(alpha**2-alpha*(1._kp+2._kp* &
                                       beta*gamma_lmi*x**gamma_lmi)+beta*gamma_lmi* &
                                       x**gamma_lmi*(1._kp+gamma_lmi*(-1._kp+beta*x**gamma_lmi)))

  end function lmi3_norm_deriv_second_potential



!epsilon_one(x)
  function lmi3_epsilon_one(x,gamma_lmi,M)    
    implicit none
    real(kp) :: lmi3_epsilon_one
    real(kp), intent(in) :: x,gamma_lmi,M

    real(kp) ::alpha,beta
    alpha=lmi3_alpha(gamma_lmi,M)
    beta=lmi3_beta(gamma_lmi,M)

    lmi3_epsilon_one = (alpha-beta*gamma_lmi*x**gamma_lmi)**2/(2._kp*x**2)
    
  end function lmi3_epsilon_one


!epsilon_two(x)
  function lmi3_epsilon_two(x,gamma_lmi,M)    
    implicit none
    real(kp) :: lmi3_epsilon_two
    real(kp), intent(in) :: x,gamma_lmi,M

    real(kp) ::alpha,beta
    alpha=lmi3_alpha(gamma_lmi,M)
    beta=lmi3_beta(gamma_lmi,M)
    
    lmi3_epsilon_two = 2._kp*(alpha+beta*(-1._kp+gamma_lmi) &
                       *gamma_lmi*x**gamma_lmi)/(x**2)
    
  end function lmi3_epsilon_two


!epsilon_three(x)
  function lmi3_epsilon_three(x,gamma_lmi,M)    
    implicit none
    real(kp) :: lmi3_epsilon_three
    real(kp), intent(in) :: x,gamma_lmi,M

    real(kp) ::alpha,beta
    alpha=lmi3_alpha(gamma_lmi,M)
    beta=lmi3_beta(gamma_lmi,M)
    
    lmi3_epsilon_three = ((alpha-beta*gamma_lmi*x**gamma_lmi)*(2._kp*alpha- & 
                          beta*(-2._kp+gamma_lmi)*(-1._kp+gamma_lmi)*gamma_lmi* &
                          x**gamma_lmi))/(x**2*(alpha+beta*(-1._kp+gamma_lmi)* & 
                          gamma_lmi*x**gamma_lmi))
    
  end function lmi3_epsilon_three



!this is integral[V(phi)/V'(phi) dphi]
  function lmi3_efold_primitive(x,gamma_lmi,M)
    implicit none
    real(kp), intent(in) :: x,gamma_lmi,M
    real(kp) :: lmi3_efold_primitive

    real(kp) ::alpha,beta
    alpha=lmi3_alpha(gamma_lmi,M)
    beta=lmi3_beta(gamma_lmi,M)

    if (alpha.eq.0._kp) stop 'lmi3_efold_primitive: alpha=0!'
    if (gamma_lmi.eq.0._kp) stop 'lmi3_efold_primitive: gamma=0!'

    lmi3_efold_primitive = x**2/(2._kp*alpha)*hypergeom_2F1(1._kp,2._kp/gamma_lmi, &
                           2._kp/(gamma_lmi+1._kp),beta*gamma_lmi/alpha*x**gamma_lmi)

    !Approximated trajectory in the x/x_0 >>1 limit
    lmi3_efold_primitive = 1._kp/(beta*gamma_lmi*(gamma_lmi-2._kp))*x**(2._kp-gamma_lmi)

  end function lmi3_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function lmi3_x_trajectory(bfold,xend,gamma_lmi,M)
    implicit none
    real(kp), intent(in) :: bfold, gamma_lmi, xend,M
    real(kp) :: lmi3_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi3Data

  
    real(kp) ::alpha,beta
    alpha=lmi3_alpha(gamma_lmi,M)
    beta=lmi3_beta(gamma_lmi,M)

    maxi = xend*(1._kp-epsilon(1._kp))
    mini = (alpha/(beta*gamma_lmi))**(1._kp/gamma_lmi)*(1._kp+epsilon(1._kp))
  

    lmi3Data%real1 = gamma_lmi
    lmi3Data%real2 = M
    lmi3Data%real3 = -bfold + lmi3_efold_primitive(xend,gamma_lmi,M)
    
    lmi3_x_trajectory = zbrent(find_lmi3traj,mini,maxi,tolFind,lmi3Data)
       
  end function lmi3_x_trajectory

  function find_lmi3traj(x,lmi3Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi3Data
    real(kp) :: find_lmi3traj
    real(kp) :: gamma_lmi,M,NplusNuend

    gamma_lmi = lmi3Data%real1
    M = lmi3Data%real2
    NplusNuend = lmi3Data%real2

    find_lmi3traj = lmi3_efold_primitive(x,gamma_lmi,M) - NplusNuend
   
  end function find_lmi3traj



end module lmi3sr
