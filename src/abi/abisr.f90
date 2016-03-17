!Exact hubble flow functions and potentials for Alpha-Beta Inflation
!
!For alpha not equal to 2 and not equal to 1:
!
!V(phi) = M^4 {1 - (1/2) beta/[(1-alpha/2)^2 x^2/(3 beta)]^[alpha/(2-alpha)]} 
! * exp{ 3 beta/(1-alpha) * [(1-alpha/2)^2 x^2 / (3 beta)]^[(1-alpha)/(2-alpha)]}
!
!For alpha=1
!
!V(phi) = M^4 (x)^(6 beta) [ 1-6 beta^2/x^2 ]
!
!For alpha=2
!
!V(phi) = M^4 {1 - (1/2) beta exp[-2 x/sqrt(3 beta)]} exp{-3 beta exp[-x/sqrt(3 beta)]}
!
!x = phi/Mpl

module abisr
  use infprec, only : kp

  implicit none

  private

  public abi_norm_potential, abi_norm_deriv_potential
  public abi_norm_deriv_second_potential
  public abi_epsilon_one, abi_epsilon_two, abi_epsilon_three
  public abi_efold_primitive, abi_x_trajectory, abi_x_endinf
  public abi_numacc_betamax

contains


  function abi_numacc_betamax(efold,alpha)
    implicit none
    real(kp) :: abi_numacc_betamax
    real(kp), intent(in) :: efold, alpha

    real(kp), parameter :: lnPotNumAccMax = log(huge(1._kp)*epsilon(1._kp))

    if (alpha.ge.1._kp) then
       abi_numacc_betamax = huge(1._kp)
       return
    endif

    abi_numacc_betamax = lnPotNumAccMax &
         * (1._kp - alpha)/3._kp/efold**(1._kp-alpha)

  end function abi_numacc_betamax




  function abi_norm_potential(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: abi_norm_potential

    real(kp) :: y
    
    if (alpha.eq.1._kp) then

       y = x*x

       abi_norm_potential = y**(3._kp*beta) &
            *(1._kp-6._kp*beta*beta/y)
       return

    elseif (alpha.eq.2._kp) then

       y = exp(-x/sqrt(3._kp*beta))

       abi_norm_potential = (1._kp - 0.5_kp*beta*y*y) &
            * exp(-3._kp*beta*y)
       return

    endif

    y = (1._kp-0.5_kp*alpha)**2/(3._kp*beta)*x*x

    abi_norm_potential = (1._kp - 0.5_kp*beta/y**(alpha/(2._kp-alpha))) &
         * exp(3._kp*beta/(1._kp-alpha)*y**((1._kp-alpha)/(2._kp-alpha)))

  end function abi_norm_potential



  function abi_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: abi_norm_deriv_potential

    real(kp) :: y

    if (alpha.eq.1_kp) then

       y = x*x

       abi_norm_deriv_potential = 6._kp*beta*y**(3._kp*beta-1.5_kp) &
            * (y + 2._kp - 6._kp*beta)
       return

    elseif (alpha.eq.2._kp) then

       y = exp(-x/sqrt(3._kp*beta))

       abi_norm_deriv_potential = sqrt(beta/12._kp)*y*exp(-3._kp*beta*y) &
            * (-3._kp*beta*y*y + 2._kp*y + 6._kp)
       return

    endif

    y = (1._kp-0.5_kp*alpha)**2/(3._kp*beta)*x*x

    abi_norm_deriv_potential = -((-2._kp + alpha)*x*y**(1._kp/(-2._kp + alpha)) &
         * (6._kp + alpha*y**(1._kp/(-2._kp + alpha))  &
         - 3._kp*beta*y**(alpha/(-2._kp + alpha)))) &
         /(12._kp*exp((3*beta*y**((-1._kp + alpha)/(-2._kp + alpha)))/(-1._kp + alpha)))

  end function abi_norm_deriv_potential




  function abi_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: abi_norm_deriv_second_potential

    real(kp) :: y

    if (alpha.eq.1._kp) then

       y = x*x

       abi_norm_deriv_second_potential = (-6._kp*beta*y**(3._kp*beta) &
            *(6._kp - 30._kp*beta + 36._kp*beta**2 + y - 6._kp*beta*y))/y/y
       return

    elseif (alpha.eq.2._kp) then

       y = exp(-x/sqrt(3._kp*beta))

       abi_norm_deriv_second_potential = -(y*(6._kp + (4._kp - 18._kp*beta)*y &
            - 15._kp*beta*y**2 + 9._kp*beta**2*y**3))/(6._kp*exp(3*beta*y))
       return

    endif

    y = (1._kp-0.5_kp*alpha)**2/(3._kp*beta)*x*x

    abi_norm_deriv_second_potential = (y**(2._kp/(-2._kp + alpha)) &
         *(-alpha**2 + alpha*(-2._kp - 6._kp*y**(1._kp/(2._kp - alpha)) &
         + 15._kp*beta*y**((-1._kp + alpha)/(-2._kp + alpha))) &
         - 18._kp*beta*(-2._kp*y + beta*y**((2*(-1._kp + alpha))/(-2._kp + alpha))))) &
         /(12._kp*exp((3._kp*beta*y**((-1._kp + alpha)/(-2._kp + alpha)))/(-1._kp + alpha)))

  end function abi_norm_deriv_second_potential




!first hubble flow function (eps1H, exact)
  function abi_epsilon_one(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: abi_epsilon_one

    real(kp) :: y

    if (alpha.eq.2._kp) then
       abi_epsilon_one = 1.5_kp*beta*exp(-2._kp*x/sqrt(3._kp*beta))
       return
    endif

    y = (1._kp-alpha/2._kp)**2/(3._kp*beta)*x*x

    abi_epsilon_one = 1.5_kp*beta/y**(alpha/(2._kp-alpha))
        
  end function abi_epsilon_one




!second hubble flow function (eps2H, exact)
  function abi_epsilon_two(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: abi_epsilon_two

    real(kp) :: y

    if (alpha.eq.2._kp) then
       abi_epsilon_two = 2._kp*exp(-x/sqrt(3._kp*beta))
       return
    endif

    y = (1._kp-alpha/2._kp)**2/(3._kp*beta)*x*x

    abi_epsilon_two = alpha/y**(1._kp/(2._kp-alpha))

  end function abi_epsilon_two
 



!third hubble flow function (eps3H, exact)
  function abi_epsilon_three(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: abi_epsilon_three

    abi_epsilon_three = abi_epsilon_two(x,alpha,beta)/alpha

  end function abi_epsilon_three



!returns the field value at which inflation ends
  function abi_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: abi_x_endinf
    
    if (alpha.eq.2._kp) then
       abi_x_endinf = 0.5_kp*sqrt(3._kp*beta)*log(1.5_kp*beta)
       return
    endif

    abi_x_endinf = (1.5_kp*beta)**(1._kp/alpha) * sqrt(2._kp)/abs(1._kp-0.5_kp*alpha)

  end function abi_x_endinf




  function abi_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: abi_efold_primitive

    real(kp) :: y, bfold

    if (alpha.eq.2._kp) then
       bfold = sqrt(1.5_kp*beta) - exp(x/sqrt(3._kp*beta))
       abi_efold_primitive = - bfold
       return
    endif

    y = (1._kp-alpha/2._kp)**2/(3._kp*beta)*x*x

    bfold = (1.5_kp*beta)**(1._kp/alpha) - y**(1._kp/(2._kp-alpha))

    abi_efold_primitive = -bfold      

  end function abi_efold_primitive




  function abi_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,beta
    real(kp) :: abi_x_trajectory

    if (alpha.eq.2._kp) then
       abi_x_trajectory = sqrt(3._kp*beta)*log(sqrt(1.5_kp*beta) - bfold)
       return
    endif

   abi_x_trajectory = sqrt(3._kp*beta)/(1._kp-0.5_kp*alpha) &
        * ((1.5_kp*beta)**(1._kp/alpha) - bfold)**(1._kp-0.5_kp*alpha)

  end function abi_x_trajectory
  


end module abisr
