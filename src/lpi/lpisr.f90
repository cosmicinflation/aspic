!slow-roll functions for the logarithmic potential inflation models
!
!V(phi) = M**4 x**p log(x)**q
!
!x = phi/phi0
!phi0 = phi0/Mp

module lpisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : ei

  implicit none

  private

  public lpi_norm_potential, lpi_epsilon_one, lpi_epsilon_two, lpi_epsilon_three
  public lpi_norm_deriv_potential, lpi_norm_deriv_second_potential
  public lpi_x_endinf, lpi_efold_primitive, lpi_x_trajectory


 
contains
!returns V/M**4
  function lpi_norm_potential(x,p,q,phi0)
    implicit none
    real(kp) :: lpi_norm_potential
    real(kp), intent(in) :: x,phi0,p,q

    lpi_norm_potential = x**p*log(x)**q

  end function lpi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function lpi_norm_deriv_potential(x,p,q,phi0)
    implicit none
    real(kp) :: lpi_norm_deriv_potential
    real(kp), intent(in) :: x,phi0,p,q

   lpi_norm_deriv_potential = x**(-1._kp+p)*log(x)**(-1._kp+q)*(q+p*log(x))

  end function lpi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function lpi_norm_deriv_second_potential(x,p,q,phi0)
    implicit none
    real(kp) :: lpi_norm_deriv_second_potential
    real(kp), intent(in) :: x,phi0,p,q

    lpi_norm_deriv_second_potential = x**(-2._kp+p)*log(x)**(-2._kp+q)*((-1._kp+q)*q+ &
                                      (-1._kp+2._kp*p)*q*log(x)+(-1._kp+p)*p*log(x)**2)

  end function lpi_norm_deriv_second_potential



!epsilon_one(x)
  function lpi_epsilon_one(x,p,q,phi0)    
    implicit none
    real(kp) :: lpi_epsilon_one
    real(kp), intent(in) :: x,phi0,p,q
    
    lpi_epsilon_one = (q+p*log(x))**2/(2._kp*phi0**2*x**2*log(x)**2)
    
  end function lpi_epsilon_one


!epsilon_two(x)
  function lpi_epsilon_two(x,p,q,phi0)    
    implicit none
    real(kp) :: lpi_epsilon_two
    real(kp), intent(in) :: x,phi0,p,q
    
    lpi_epsilon_two = (2._kp*(q+log(x)*(q+p*log(x))))/(phi0**2*x**2*log(x)**2)
    
  end function lpi_epsilon_two


!epsilon_three(x)
  function lpi_epsilon_three(x,p,q,phi0)    
    implicit none
    real(kp) :: lpi_epsilon_three
    real(kp), intent(in) :: x,phi0,p,q
    
    lpi_epsilon_three = ((q+p*log(x))*(2._kp*q+log(x)*(3._kp*q+2._kp*log(x)*(q+p*log(x)))))/ &
                        (phi0**2*x**2*log(x)**2*(q+log(x)*(q+p*log(x))))

    
  end function lpi_epsilon_three


!returns the value for x=phi/Mp defined as epsilon1=1, where inflation ends
  function lpi_x_endinf(p,q,phi0)
    implicit none
    real(kp), intent(in) :: phi0,p,q
    real(kp) :: lpi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lpiData

    mini =1._kp+q/sqrt(2._kp)*sqrt(epsilon(1._kp))/phi0 ! uses an asymptotic expression for epsilon1 when x->1 and requiring eps1<1/epsilon(1._kp) to avoid numerical issues
    maxi = max(p*10._kp**(3)/(sqrt(2._kp)*phi0),10._kp**(3)) ! uses an asymptotic expression for epsilon1 when x->Infinity and and set the upper bound such that eps1 ~ 10^(-6) there

    lpiData%real1 = p
    lpiData%real2 = q
    lpiData%real3 = phi0
    
    lpi_x_endinf = zbrent(find_lpi_x_endinf,mini,maxi,tolFind,lpiData)
   
   
  end function lpi_x_endinf

  function find_lpi_x_endinf(x,lpiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lpiData
    real(kp) :: find_lpi_x_endinf
    real(kp) :: phi0,p,q

    p = lpiData%real1
    q = lpiData%real2
    phi0 = lpiData%real3

    find_lpi_x_endinf = lpi_epsilon_one(x,p,q,phi0)-1._kp
   
  end function find_lpi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function lpi_efold_primitive(x,p,q,phi0)
    implicit none
    real(kp), intent(in) :: x,phi0,p,q
    real(kp) :: lpi_efold_primitive

    lpi_efold_primitive = phi0**2*(x**2/(2._kp*p)- &
                          q/(p**2)*exp(-2._kp*q/p)*ei(2._kp*q/p+2._kp*log(x)))

  end function lpi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function lpi_x_trajectory(bfold,xend,p,q,phi0)
    implicit none
    real(kp), intent(in) :: bfold, phi0,p,q, xend
    real(kp) :: lpi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lpiData

  
    mini = xend*(1._kp+epsilon(1._kp))
    maxi = sqrt(xend**2+2._kp*p/(phi0**2)*10._kp**(4.)) !Uses an approximate solution for the slow roll trajectory when x>>1 and sets maxi at ~ 10^4 efolds away from xend

    lpiData%real1 = p
    lpiData%real2 = q
    lpiData%real3 = phi0
    lpiData%real4 = -bfold + lpi_efold_primitive(xend,p,q,phi0)
    
    lpi_x_trajectory = zbrent(find_lpitraj,mini,maxi,tolFind,lpiData)
       
  end function lpi_x_trajectory

  function find_lpitraj(x,lpiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lpiData
    real(kp) :: find_lpitraj
    real(kp) :: phi0,p,q,NplusNuend

    p= lpiData%real1
    q= lpiData%real2
    phi0= lpiData%real3
    NplusNuend = lpiData%real4

    find_lpitraj = lpi_efold_primitive(x,p,q,phi0) - NplusNuend
   
  end function find_lpitraj

  
end module lpisr
