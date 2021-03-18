!slow-roll functions for the double exponential inflation potential
!
!V(phi) = M**4 * [ exp(beta x) - beta**2 exp(x/beta) ]
!
!x = phi/phi0
!beta is a free positive parameter

module deisr
  use infprec, only : kp,tolkp,toldp,transfert
  use inftools, only : zbrent, easydverk
  implicit none

  private

  public  dei_norm_potential, dei_epsilon_one, dei_epsilon_two, dei_epsilon_three
  public  dei_x_endinf, dei_efold_primitive, dei_x_trajectory
  public  dei_norm_deriv_potential, dei_norm_deriv_second_potential
 
contains
!returns V/M**4
  function dei_norm_potential(x,beta,phi0)
    implicit none
    real(kp) :: dei_norm_potential
    real(kp), intent(in) :: x,beta
    real(kp), intent(in) :: phi0

    dei_norm_potential = exp(beta*x)-beta**2*exp(x/beta)

  end function dei_norm_potential


!returns the first derivative of the potential with respect to x=phi/phi0, divided by M**4
  function dei_norm_deriv_potential(x,beta,phi0)
    implicit none
    real(kp) :: dei_norm_deriv_potential
    real(kp), intent(in) :: x,beta
    real(kp), intent(in) :: phi0

   dei_norm_deriv_potential = beta*exp(beta*x)-beta*exp(x/beta)

  end function dei_norm_deriv_potential



!returns the second derivative of the potential with respect to x=phi/phi0, divided by M**4
  function dei_norm_deriv_second_potential(x,beta,phi0)
    implicit none
    real(kp) :: dei_norm_deriv_second_potential
    real(kp), intent(in) :: x,beta
    real(kp), intent(in) :: phi0

    dei_norm_deriv_second_potential = beta**2*exp(beta*x)-exp(x/beta)

  end function dei_norm_deriv_second_potential



!epsilon_one(x)
  function dei_epsilon_one(x,beta,phi0)    
    implicit none
    real(kp) :: dei_epsilon_one
    real(kp), intent(in) :: x,beta,phi0

  
    dei_epsilon_one = (beta**2*(exp(x/beta)-exp(beta*x))**2)/(2._kp*(-(beta**2*exp(x/beta)) &
		      +exp(beta*x))**2*phi0**2)

    
  end function dei_epsilon_one


!epsilon_two(x)
  function dei_epsilon_two(x,beta,phi0)    
    implicit none
    real(kp) :: dei_epsilon_two
    real(kp), intent(in) :: x,beta,phi0
    
    dei_epsilon_two = (2*(-1._kp+beta**2)**2*exp(x/beta+beta*x))/((-(beta**2*exp(x/beta)) &
		      +exp(beta*x))**2*phi0**2)
   
  end function dei_epsilon_two


!epsilon_three(x)
  function dei_epsilon_three(x,beta,phi0)    
    implicit none
    real(kp) :: dei_epsilon_three
    real(kp), intent(in) :: x,beta,phi0
    
    dei_epsilon_three = ((-1._kp+ beta**2)*(-exp(x/beta)+exp(beta*x))*(beta**2*exp(x/beta) &
			+exp(beta*x)))/((-(beta**2*exp(x/beta))+exp(beta*x))**2*phi0**2)
    
  end function dei_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function dei_x_endinf(beta,phi0)
    implicit none
    real(kp), intent(in) :: beta,phi0
    real(kp) :: dei_x_endinf
    
    dei_x_endinf = (beta*log((beta*(1._kp+sqrt(2._kp)*beta*phi0))/(beta+sqrt(2._kp)*phi0)))/(-1._kp+beta**2)

  end function dei_x_endinf

!this is integral[V(phi)/V'(phi) dphi]
  function dei_efold_primitive(x,beta,phi0)
    implicit none
    real(kp), intent(in) :: x,beta,phi0
    real(kp) :: dei_efold_primitive

    dei_efold_primitive = ((1._kp+beta**2)/beta*x-log(exp(x/beta)-exp(beta*x)))*phi0**2

  end function dei_efold_primitive

  subroutine find_dei_efold_primitive(n,x,y,yprime,deiData)
    implicit none          
    integer :: n
    real(kp) :: x
    real(kp), dimension(n) :: y, yprime
    type(transfert), optional, intent(inout) :: deiData
    real(kp) :: beta, phi0, x4p2n, phi4p2n

    beta = deiData%real1
    phi0 = deiData%real2

    yprime(1) = phi0**2*((1.-x**2)**2+beta*x**4*(log(x)-0.25_kp)+ &
                beta/4._kp)/(4._kp*x*(-1._kp+x**2+x**2*beta*log(x)))

  end subroutine find_dei_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function dei_x_trajectory(bfold,xend,beta,phi0)
    implicit none
    real(kp), intent(in) :: bfold, beta, phi0, xend
    real(kp) :: dei_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: deiData

  
    mini = tolkp
    maxi = dei_x_endinf(beta,phi0)*(1._kp-epsilon(1._kp))
  

    deiData%real1 = beta
    deiData%real2 = phi0	
    deiData%real3 = -bfold + dei_efold_primitive(xend,beta,phi0)
    
    dei_x_trajectory = zbrent(find_dei_x_trajectory,mini,maxi,tolFind,deiData)
       
  end function dei_x_trajectory

  function find_dei_x_trajectory(x,deiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: deiData
    real(kp) :: find_dei_x_trajectory
    real(kp) :: beta,phi0,NplusNuend

    beta = deiData%real1
    phi0 = deiData%real2
    NplusNuend = deiData%real3

    find_dei_x_trajectory = dei_efold_primitive(x,beta,phi0) - NplusNuend
   
  end function find_dei_x_trajectory



end module deisr
