!slow-roll functions for the pure arctan inflation potential
!
!V(phi) = M^4 arctan(x)
!
!x = phi/mu

module paisr
  use infprec, only : pi, kp, tolkp, transfert
  use inftools, only : zbrent
  implicit none

  private

  public pai_norm_potential, pai_epsilon_one, pai_epsilon_two, pai_epsilon_three
  public pai_efold_primitive, pai_x_trajectory
  public pai_norm_deriv_potential, pai_norm_deriv_second_potential
  public pai_numacc_xinimax, pai_x_endinf
  
contains


!returns V/M^4 as a function of x=phi/mu
  function pai_norm_potential(x,mu)
    implicit none
    real(kp) :: pai_norm_potential
    real(kp), intent(in) :: x,mu

    pai_norm_potential = atan(x)

  end function pai_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function pai_norm_deriv_potential(x,mu)
    implicit none
    real(kp) :: pai_norm_deriv_potential
    real(kp), intent(in) :: x,mu

   pai_norm_deriv_potential = 1._kp/(1._kp + x*x)

  end function pai_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function pai_norm_deriv_second_potential(x,mu)
    implicit none
    real(kp) :: pai_norm_deriv_second_potential
    real(kp), intent(in) :: x,mu

    pai_norm_deriv_second_potential = -2._kp*x/(1._kp + x*x)**2

  end function pai_norm_deriv_second_potential


!how far we can start inflation and having epsilon1 > machine precision  
  function pai_numacc_xinimax(mu)
    implicit none
    real(kp) :: pai_numacc_xinimax
    real(kp), intent(in) :: mu

    real(kp), parameter :: sqrsqrepsmin = epsilon(1._kp)**(0.25_kp)
    
    pai_numacc_xinimax = 1._kp/(sqrt(pi*mu)*sqrsqrepsmin)
    
  end function pai_numacc_xinimax
  


!epsilon_one(x)
  function pai_epsilon_one(x,mu)    
    implicit none
    real(kp) :: pai_epsilon_one
    real(kp), intent(in) :: x,mu
    
    pai_epsilon_one = 0.5_kp/(atan(x)**2*(1._kp+x*x)**2)/mu/mu
    
  end function pai_epsilon_one


!epsilon_two(x)
  function pai_epsilon_two(x,mu)    
    implicit none
    real(kp) :: pai_epsilon_two
    real(kp), intent(in) :: x,mu
    
    pai_epsilon_two = (2._kp + 4._kp*x*atan(x)) &
         /(mu**2 * (1 + x*x)**2 * atan(x)**2)

  end function pai_epsilon_two


!epsilon_three(x)
  function pai_epsilon_three(x,mu)    
    implicit none
    real(kp) :: pai_epsilon_three
    real(kp), intent(in) :: x,mu
    
    pai_epsilon_three = (2._kp + 6._kp*x*atan(x) &
         + (-2._kp + 6._kp*x**2)*atan(x)**2) &
         /(mu**2 * (1._kp + x*x)**2*atan(x)**2 &
         * (1._kp + 2._kp*x*atan(x)))
    
  end function pai_epsilon_three

 


!returns x at the end of inflation defined as epsilon1=1
  function pai_x_endinf(mu)
    implicit none
    real(kp), intent(in) :: mu
    real(kp) :: pai_x_endinf
    real(kp) :: mini, maxi
    real(kp), parameter :: tolFind=tolkp
    type(transfert) :: paiData

    mini = 0._kp
    maxi = pai_numacc_xinimax(mu)

    paiData%real1 = mu

    pai_x_endinf = zbrent(find_pai_x_endinf,mini,maxi,tolFind,paiData)
    
  end function pai_x_endinf

  function find_pai_x_endinf(x,paiData)
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: paiData
    real(kp) :: find_pai_x_endinf
    real(kp) :: mu

    mu = paiData%real1

    find_pai_x_endinf = mu*sqrt(2._kp+2._kp*x*x)*atan(x) - 1._kp
    
  end function find_pai_x_endinf

    

!this is integral[V(phi)/V'(phi) dphi]
  function pai_efold_primitive(x,mu)
    implicit none
    real(kp), intent(in) :: x,mu
    real(kp) :: pai_efold_primitive


    pai_efold_primitive = (mu**2 * (-x**2 + 2._kp*x*(3._kp + x**2)*atan(x) &
         - 2._kp*log(1._kp + x**2)))/6._kp
    
  end function pai_efold_primitive


  
!returns x at bfold=-efolds before the end of inflation
  function pai_x_trajectory(bfold,xend,mu)
    implicit none
    real(kp), intent(in) :: bfold,xend,mu
    real(kp) :: pai_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: paiData

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = pai_numacc_xinimax(mu)

    paiData%real1 = mu
    paiData%real2 = -bfold + pai_efold_primitive(xend,mu)

    pai_x_trajectory = zbrent(find_pai_x_trajectory,mini,maxi,tolFind,paiData)

  end function pai_x_trajectory

  function find_pai_x_trajectory(x,paiData)
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: paiData
    real(kp) :: find_pai_x_trajectory
    real(kp) :: mu,NplusPrimEnd

    mu = paiData%real1
    NplusPrimEnd = paiData%real2

    find_pai_x_trajectory = pai_efold_primitive(x,mu) - NplusPrimEnd

  end function find_pai_x_trajectory



  
end module paisr
