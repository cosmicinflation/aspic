!slow-roll functions for the logarithmic potential inflation models
!
!V(phi) = M**4 x**p log(x)**q
!
!x = phi/phi0
!phi0 = phi0/Mp

module lpicommon
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : ei

  implicit none
  
  private

  real(kp), parameter :: xBig = 1._kp/(10._kp*epsilon(1._kp))
  
  public xBig
  public lpi_norm_potential, lpi_epsilon_one, lpi_epsilon_two, lpi_epsilon_three
  public lpi_norm_deriv_potential, lpi_norm_deriv_second_potential
  public lpi_x_potmax, lpi_efold_primitive, lpi_x_epstwozero
  public lpi23_sanity_check, find_lpi_x_trajectory

 
contains

!q should be integer and even
  subroutine lpi23_sanity_check(p,q,phi0)
    implicit none
    real(kp), intent(in) :: q
    real(kp), intent(in), optional :: p,phi0
    real(kp), parameter :: epstwominMax = 10._kp

    real(kp) :: vevMin

    if (modulo(q,2._kp).ne.0._kp) then
       write(*,*)'q= ',q
       stop 'lpi23_sanity_check: q should be positive and even!'
    endif
    
    if (present(p).and.present(phi0)) then
!this is phi0 such that the min of eps2 in the domain lpi23,
!i.e. eps2(xVmax)=epstwominMax
       vevMin = p * exp(q/p) * sqrt(2._kp/q/epstwominMax)

       if (phi0.lt.vevMin) then
          write(*,*)'vevmin= phi0= ',vevMin,phi0
          stop 'lpi23_sanity_check: min(eps2) >> 1, phi0 too small!'
       end if

    end if

  end subroutine lpi23_sanity_check




!returns V/M**4
  function lpi_norm_potential(x,p,q,phi0)
    implicit none
    real(kp) :: lpi_norm_potential
    real(kp), intent(in) :: x,p,q
    real(kp), intent(in), optional :: phi0

    lpi_norm_potential = x**p*log(x)**q

  end function lpi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function lpi_norm_deriv_potential(x,p,q,phi0)
    implicit none
    real(kp) :: lpi_norm_deriv_potential
    real(kp), intent(in) :: x,p,q
    real(kp), intent(in), optional :: phi0

   lpi_norm_deriv_potential = x**(-1._kp+p)*log(x)**(-1._kp+q)*(q+p*log(x))

  end function lpi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function lpi_norm_deriv_second_potential(x,p,q,phi0)
    implicit none
    real(kp) :: lpi_norm_deriv_second_potential
    real(kp), intent(in) :: x,p,q
    real(kp), intent(in), optional :: phi0

    lpi_norm_deriv_second_potential = x**(-2._kp+p)*log(x)**(-2._kp+q)*((-1._kp+q)*q+ &
         (-1._kp+2._kp*p)*q*log(x)+(-1._kp+p)*p*log(x)**2)

  end function lpi_norm_deriv_second_potential



!epsilon_one(x)
  function lpi_epsilon_one(x,p,q,phi0)    
    implicit none
    real(kp) :: lpi_epsilon_one
    real(kp), intent(in) :: x,p,q
    real(kp), intent(in) :: phi0

    lpi_epsilon_one = (q+p*log(x))**2/(2._kp*phi0**2*x**2*log(x)**2)
    
  end function lpi_epsilon_one


!epsilon_two(x)
  function lpi_epsilon_two(x,p,q,phi0)    
    implicit none
    real(kp) :: lpi_epsilon_two
    real(kp), intent(in) :: x,p,q
    real(kp), intent(in) :: phi0

    lpi_epsilon_two = (2._kp*(q+log(x)*(q+p*log(x))))/(phi0**2*x**2*log(x)**2)
    
  end function lpi_epsilon_two


!epsilon_three(x)
  function lpi_epsilon_three(x,p,q,phi0)    
    implicit none
    real(kp) :: lpi_epsilon_three
    real(kp), intent(in) :: x,p,q
    real(kp), intent(in) :: phi0

    lpi_epsilon_three = ((q+p*log(x))*(2._kp*q+log(x)*(3._kp*q+2._kp*log(x)*(q+p*log(x)))))/ &
         (phi0**2*x**2*log(x)**2*(q+log(x)*(q+p*log(x))))

    
  end function lpi_epsilon_three


  function lpi_x_potmax(p,q,phi0)
    implicit none
    real(kp), intent(in) :: p,q
    real(kp), intent(in) :: phi0
    real(kp) :: lpi_x_potmax

    lpi_x_potmax = exp(-q/p)
    
  end function lpi_x_potmax

  

  function lpi_x_epstwozero(p,q,phi0)
    implicit none
    real(kp), intent(in) :: p,q
    real(kp), intent(in) :: phi0
    real(kp), dimension(2) :: lpi_x_epstwozero

    lpi_x_epstwozero(1) = exp( (-q - sqrt(-4._kp*p*q + q*q))/(2._kp*p) )
    lpi_x_epstwozero(2) = exp( (-q + sqrt(-4._kp*p*q + q*q))/(2._kp*p) )
        
  end function lpi_x_epstwozero



  function lpi_epstwo_potmax(p,q,phi0)
    implicit none
    real(kp), intent(in) :: p,q, phi0
    real(kp) :: lpi_epstwo_potmax
    real(kp) :: xVmax

    xVmax = lpi_x_potmax(p,q,phi0)
    
    lpi_epstwo_potmax = lpi_epsilon_two(xVmax,p,q,phi0)

  end function lpi_epstwo_potmax


!this is integral[V(phi)/V'(phi) dphi]
  function lpi_efold_primitive(x,p,q,phi0)
    implicit none
    real(kp), intent(in) :: x,phi0,p,q
    real(kp) :: lpi_efold_primitive

    if (x.le.0._kp) stop 'lpi_efold_primitive: x <=0'
    
    lpi_efold_primitive = phi0**2*(x**2/(2._kp*p)- &
         q/(p**2)*exp(-2._kp*q/p)*ei(2._kp*q/p+2._kp*log(x)))

  end function lpi_efold_primitive




  function find_lpi_x_trajectory(x,lpiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lpiData
    real(kp) :: find_lpi_x_trajectory
    real(kp) :: phi0,p,q,NplusNuend

    p= lpiData%real1
    q= lpiData%real2
    phi0= lpiData%real3
    NplusNuend = lpiData%real4

    find_lpi_x_trajectory = lpi_efold_primitive(x,p,q,phi0) - NplusNuend
   
  end function find_lpi_x_trajectory

  
end module lpicommon
