!slow-roll functions for the S-dual inflation potential
!
!V(phi) = M^4 sech(x)
!
!x = phi/phi0

module sdisr
  use infprec, only : kp
  use cosmopar, only : efoldNum
  implicit none

  private

  public  sdi_norm_potential, sdi_epsilon_one, sdi_epsilon_two, sdi_epsilon_three
  public  sdi_efold_primitive, sdi_x_trajectory, sdi_check_params
  public  sdi_norm_deriv_potential, sdi_norm_deriv_second_potential
 
contains

! Hardprior condition on phi0 and xend
  function sdi_check_params(phi0, xend)
    implicit none
    logical :: sdi_check_params
    real(kp), intent(in) :: phi0, xend
    real(kp) :: Nmax, phi0min, xendmax

    Nmax = sdi_efold_primitive(epsilon(1._kp),phi0)-sdi_efold_primitive(xend,phi0)    ! Virtually the max number of efolds is infinite, but this is the maximum number of efolds one can numerically resolve, given the accuracy (kp)
    xendmax = -log(epsilon(1._kp)) ! this the maximum value the exponential of which can be dealt with given the accuracy (kp)
    phi0min = sqrt(2.)/2._kp
    phi0min = 2._kp !ruled out otherwise anyway

    sdi_check_params = ((phi0 .ge. phi0min) .and. (Nmax .ge. efoldNum) .and. (xend<xendmax) )

  end function sdi_check_params


!returns V/M^4 as a function of x=phi/phi0
  function sdi_norm_potential(x)
    implicit none
    real(kp) :: sdi_norm_potential
    real(kp), intent(in) :: x

    sdi_norm_potential = 1._kp/cosh(x)

  end function sdi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function sdi_norm_deriv_potential(x)
    implicit none
    real(kp) :: sdi_norm_deriv_potential
    real(kp), intent(in) :: x

   sdi_norm_deriv_potential = -tanh(x)/cosh(x)

  end function sdi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function sdi_norm_deriv_second_potential(x)
    implicit none
    real(kp) :: sdi_norm_deriv_second_potential
    real(kp), intent(in) :: x

    sdi_norm_deriv_second_potential = 1._kp/cosh(x)-2._kp/cosh(x)**3

  end function sdi_norm_deriv_second_potential



!epsilon_one(x)
  function sdi_epsilon_one(x,phi0)    
    implicit none
    real(kp) :: sdi_epsilon_one
    real(kp), intent(in) :: x,phi0
    
    sdi_epsilon_one = tanh(x)**2/(2._kp*phi0**2)
    
  end function sdi_epsilon_one


!epsilon_two(x)
  function sdi_epsilon_two(x,phi0)    
    implicit none
    real(kp) :: sdi_epsilon_two
    real(kp), intent(in) :: x,phi0
    
    sdi_epsilon_two = (2._kp/cosh(x)**2)/phi0**2

  end function sdi_epsilon_two


!epsilon_three(x)
  function sdi_epsilon_three(x,phi0)    
    implicit none
    real(kp) :: sdi_epsilon_three
    real(kp), intent(in) :: x,phi0
    
    sdi_epsilon_three = (-2._kp*tanh(x)**2)/phi0**2
    
  end function sdi_epsilon_three

!returns x at the end of inflation
  function sdi_x_endinf(phi0,xend)
    implicit none
    real(kp), intent(in) :: phi0,xend
    real(kp) :: sdi_x_endinf
   
    sdi_x_endinf = xend
   
  end function sdi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function sdi_efold_primitive(x,phi0)
    implicit none
    real(kp), intent(in) :: x,phi0
    real(kp) :: sdi_efold_primitive

    if (phi0.eq.0._kp) stop 'sdi_efold_primitive: phi0=0!'

    sdi_efold_primitive = -phi0**2*log(sinh(x))

  end function sdi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function sdi_x_trajectory(bfold,xend,phi0)
    implicit none
    real(kp), intent(in) :: bfold, phi0, xend
    real(kp) :: sdi_x_trajectory
    
    sdi_x_trajectory = asinh(sinh(xend)*exp(bfold/phi0**2))
       
  end function sdi_x_trajectory


  
end module sdisr
