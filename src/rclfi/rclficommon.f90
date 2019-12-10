!common functions for the radiatively corrected large field inflation potential
!
!V(phi) = M^4 [x^p + alpha x^4 ln(x) ]
!
!x = phi/mu
!
!with no assumption on alpha, but p > 0 and p <> 4
!
module rclficommon
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  use inftools, only : zbrent
  implicit none

  private

  public rclfi_check_params, rclfi_alpha_one, rclfi_alpha_zero
  public rclfi_norm_potential, rclfi_epsilon_one, rclfi_epsilon_two, rclfi_epsilon_three
  public rclfi_norm_deriv_potential, rclfi_norm_deriv_second_potential
  public rclfi_x_derivpotzero, rclfi_x_potzero
  public find_rclfi_x_endinf, find_rclfi_efold_primitive
 
contains


  
  function rclfi_check_params(alpha,p,mu)
    implicit none
    logical :: rclfi_check_params
    real(kp), intent(in) :: alpha,p,mu

    rclfi_check_params = (p.gt.0._kp).and.(p.ne.4._kp)

  end function rclfi_check_params
  
!returns the value of alpha at which the potential develops two
!extrema
  function rclfi_alpha_one(p)
    implicit none
    real(kp) :: rclfi_alpha_one
    real(kp), intent(in) :: p

    rclfi_alpha_one = -0.25_kp*p*(p-4._kp)*exp(2._kp-0.25_kp*p)

  end function rclfi_alpha_one


!returns the value of alpha at which the potential develops two
!non-vanishing zeros
  function rclfi_alpha_zero(p)
    implicit none
    real(kp) :: rclfi_alpha_zero
    real(kp), intent(in) :: p

    rclfi_alpha_zero = -(p-4._kp)*exp(1._kp)

  end function rclfi_alpha_zero


  
  !returns V/M^4
  function rclfi_norm_potential(x,alpha,p,mu)
    implicit none
    real(kp) :: rclfi_norm_potential
    real(kp), intent(in) :: x,alpha,p,mu

    rclfi_norm_potential = x**p + alpha*x**4*log(x)

  end function rclfi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function rclfi_norm_deriv_potential(x,alpha,p,mu)
    implicit none
    real(kp) :: rclfi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,p,mu

   rclfi_norm_deriv_potential = p*x**(p-1._kp) + alpha*x**3*(1._kp + 4._kp*log(x))

  end function rclfi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function rclfi_norm_deriv_second_potential(x,alpha,p,mu)
    implicit none
    real(kp) :: rclfi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,p,mu

    rclfi_norm_deriv_second_potential = p*(p-1._kp)*x**(p-2._kp) &
         + alpha*x*x*(7._kp + 12._kp*log(x))

  end function rclfi_norm_deriv_second_potential



!epsilon_one(x)
  function rclfi_epsilon_one(x,alpha,p,mu)    
    implicit none
    real(kp) :: rclfi_epsilon_one
    real(kp), intent(in) :: x,alpha,p,mu

    real(kp) :: lnx

    lnx = log(x)
    
    rclfi_epsilon_one =  (p*x**p + x**4*alpha*(1._kp + 4._kp*lnx))**2 &
         /(2._kp*x**2*mu**2*(x**p + x**4*alpha*lnx)**2)
    
  end function rclfi_epsilon_one


!epsilon_two(x)
  function rclfi_epsilon_two(x,alpha,p,mu)    
    implicit none
    real(kp) :: rclfi_epsilon_two
    real(kp), intent(in) :: x,alpha,p,mu

    real(kp) :: lnx

    lnx = log(x)
    
    rclfi_epsilon_two = (2._kp*p*x**(2._kp*p) - 2._kp*alpha*x**(4._kp + p) &
         *(7._kp - 2._kp*p +  (12._kp + (-9._kp + p)*p)*lnx) + 2._kp*alpha**2*x**8 &
         *(1._kp + lnx + 4._kp*lnx**2))/(mu**2*x**2*(x**p + alpha*x**4*lnx)**2)


  end function rclfi_epsilon_two


!epsilon_three(x)
  function rclfi_epsilon_three(x,alpha,p,mu)    
    implicit none
    real(kp) :: rclfi_epsilon_three
    real(kp), intent(in) :: x,alpha,p,mu

    real(kp) :: lnx

    lnx = log(x)
    
    rclfi_epsilon_three =  ((alpha*x**4 + p*x**p + 4._kp*alpha*x**4*lnx) &
         *(3._kp*alpha*p**2*x**(4._kp + 2._kp*p) + alpha*x**4 &
         *(2._kp*alpha**2*x**8 + 26._kp*x**(2._kp*p) - 21._kp*alpha*x**(4._kp + p)) &
         + 2._kp*p*x**p*(3._kp*alpha**2*x**8 + x**(2._kp*p) - 9._kp*alpha*x**(4._kp + p)) &
         + alpha*x**4*(3._kp*alpha**2*x**8 - (-24._kp + 20._kp*p - 9._kp*p**2 + p**3) &
         *x**(2._kp*p) - alpha*(68._kp - 30._kp*p + 3._kp*p**2) &
         *x**(4._kp + p))*lnx + alpha**2*x**8*(2._kp*alpha*x**4 &
         + (-96._kp + 74._kp*p - 15._kp*p**2 + p**3)*x**p)*lnx**2 &
         + 8._kp*alpha**3*x**12*lnx**3)) / (mu**2*x**2*(x**p + alpha*x**4*lnx)**2 &
         *(alpha**2*x**8 - 7._kp*alpha*x**(4._kp + p) + p*x**p*(2._kp*alpha*x**4 + x**p) &
         + alpha*x**4*(alpha*x**4 - (12._kp - 9._kp*p + p**2)*x**p)*lnx &
         + 4*alpha**2*x**8*lnx**2))
    
  end function rclfi_epsilon_three



!non vanishing x at which the potential is extremal
  function rclfi_x_derivpotzero(alpha,p,mu)
    implicit none
    real(kp), dimension(2) :: rclfi_x_derivpotzero
    real(kp), intent(in) :: alpha,p,mu

    real(kp) :: arg, ppm4oa,e1mpo4, xb0, xb1

    ppm4oa = 0.25_kp*p*(p-4._kp)/alpha
    e1mpo4 = exp(1._kp - 0.25_kp*p)    
    
    arg = ppm4oa*e1mpo4
    
    if (arg.eq.0._kp) stop 'rclfi_x_derivpotzero: p=4!'
    
    if (arg.gt.0._kp) then
       rclfi_x_derivpotzero = (lambert(arg,0)/ppm4oa)**(1._kp/(p-4._kp))
    elseif (arg.gt.-1._kp/exp(1._kp)) then
       xb1 = (lambert(arg,-1)/ppm4oa)**(1._kp/(p-4._kp))
       xb0 = (lambert(arg,0)/ppm4oa)**(1._kp/(p-4._kp))

       rclfi_x_derivpotzero(1) = min(xb0,xb1)
       rclfi_x_derivpotzero(2) = max(xb0,xb1)

    elseif (arg.eq.-1._kp/exp(1._kp)) then

       rclfi_x_derivpotzero = exp((8._kp-p)/(4._kp*(p-4._kp)))
       
    else
       stop 'rclfi_x_derivpotzero: [p(p-4)e^(1-p/4)]/alpha < -1/e!'
    endif

    
  end function rclfi_x_derivpotzero


!non vanishing x at which the potential vanishes
  function rclfi_x_potzero(alpha,p,mu)
    implicit none
    real(kp), dimension(2) :: rclfi_x_potzero
    real(kp), intent(in) :: alpha,p,mu
    real(kp) :: arg

    real(kp) :: xb0,xb1
    
    arg = (p-4._kp)/alpha

    if (arg.eq.0._kp) stop 'rclfi_x_potzero: p=4!'
    
    if (arg.gt.0._kp) then
       rclfi_x_potzero = (lambert(arg,0)/arg)**(1._kp/(p-4._kp))
    elseif (arg.gt.-1._kp/exp(1._kp)) then
       xb0 = (lambert(arg,0)/arg)**(1._kp/(p-4._kp))
       xb1 = (lambert(arg,-1)/arg)**(1._kp/(p-4._kp))

       rclfi_x_potzero(1) = min(xb0,xb1)
       rclfi_x_potzero(2) = max(xb0,xb1)

    elseif (arg.eq.-1._kp/exp(1._kp)) then

       rclfi_x_potzero = exp(1._kp/(p-4._kp))
       
    else
       write(*,*)'alpha= p= ',alpha,p,arg
       stop 'rclfi_x_potzero: (p-4)/alpha < -1/e!'
    endif
    
  end function rclfi_x_potzero
  

 

  function find_rclfi_x_endinf(x,rclfiData)
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: rclfiData
    real(kp) :: find_rclfi_x_endinf
    real(kp) :: alpha,p,mu
    
    alpha = rclfiData%real1
    p=rclfiData%real2
    mu=rclfiData%real3
    
    find_rclfi_x_endinf = rclfi_epsilon_one(x,alpha,p,mu) - 1._kp

  end function find_rclfi_x_endinf


  


  subroutine find_rclfi_efold_primitive(n,x,y,yprime,rclfiData)
    implicit none          
    integer :: n
    real(kp) :: x
    real(kp), dimension(n) :: y, yprime
    type(transfert), optional, intent(inout) :: rclfiData
!at order 3
    real(kp), parameter :: xtaylor = 10._kp*epsilon(1._kp)**(1._kp/3._kp)
    
    real(kp) :: xp,x4,lnx,alpha,p
    
    alpha = rclfiData%real1
    p = rclfiData%real2

    if (abs(x).lt.xtaylor) then

       yprime(1) = x/p

    else
              
       xp = x**p
       x4 = x**4
       lnx = log(x)

       yprime(1) = x*(xp+alpha*x4*lnx)/( p*xp+alpha*x4*(1._kp+4._kp*lnx) )
      
    endif
    
  end subroutine find_rclfi_efold_primitive

  

  
end module rclficommon
