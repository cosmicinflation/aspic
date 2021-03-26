!slow-roll functions for the S-dual inflation potential
!
!V(phi) = M^4 sech(x)
!
!x = phi/mu

module sdisr
  use infprec, only : kp
  implicit none

  private

  public sdi_norm_potential, sdi_epsilon_one, sdi_epsilon_two, sdi_epsilon_three
  public sdi_efold_primitive, sdi_x_trajectory
  public sdi_norm_deriv_potential, sdi_norm_deriv_second_potential
  public sdi_phizeromin, sdi_numacc_xendmin, sdi_numacc_xendmax, sdi_numacc_xinimin

contains


!to ensure eps1 < 1 over the whole domain, otherwise ruled out by eps2
!too big
  function sdi_phizeromin()
    implicit none
    real(kp) :: sdi_phizeromin

    sdi_phizeromin = sqrt(2._kp)/2._kp

  end function sdi_phizeromin


  function sdi_numacc_xinimin(mu)
    implicit none
    real(kp) :: sdi_numacc_xinimin
    real(kp), intent(in) :: mu

    sdi_numacc_xinimin = epsilon(1._kp)

  end function sdi_numacc_xinimin
    

!the minimal value of xend ensuring efoldNum at the top of the
!potential, numerical limitation only
  function sdi_numacc_xendmin(efoldNum,mu)
    implicit none
    real(kp) :: sdi_numacc_xendmin
    real(kp), intent(in) :: efoldNum, mu
    real(kp) :: xinimin

    xinimin = sdi_numacc_xinimin(mu)

    sdi_numacc_xendmin = sdi_x_trajectory(efoldNum,xinimin,mu)
    
  end function sdi_numacc_xendmin



!otherwise the potential overflows
  function sdi_numacc_xendmax(mu)
    implicit none
    real(kp) :: sdi_numacc_xendmax
    real(kp), intent(in) :: mu

    sdi_numacc_xendmax = -log(epsilon(1._kp)) 

  end function sdi_numacc_xendmax



!returns V/M^4 as a function of x=phi/mu
  function sdi_norm_potential(x,mu)
    implicit none
    real(kp) :: sdi_norm_potential
    real(kp), intent(in) :: x,mu

    sdi_norm_potential = 1._kp/cosh(x)

  end function sdi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function sdi_norm_deriv_potential(x,mu)
    implicit none
    real(kp) :: sdi_norm_deriv_potential
    real(kp), intent(in) :: x,mu

   sdi_norm_deriv_potential = -tanh(x)/cosh(x)

  end function sdi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function sdi_norm_deriv_second_potential(x,mu)
    implicit none
    real(kp) :: sdi_norm_deriv_second_potential
    real(kp), intent(in) :: x,mu

    sdi_norm_deriv_second_potential = 1._kp/cosh(x)-2._kp/cosh(x)**3

  end function sdi_norm_deriv_second_potential



!epsilon_one(x)
  function sdi_epsilon_one(x,mu)    
    implicit none
    real(kp) :: sdi_epsilon_one
    real(kp), intent(in) :: x,mu
    
    sdi_epsilon_one = tanh(x)**2/(2._kp*mu**2)
    
  end function sdi_epsilon_one


!epsilon_two(x)
  function sdi_epsilon_two(x,mu)    
    implicit none
    real(kp) :: sdi_epsilon_two
    real(kp), intent(in) :: x,mu
    
    sdi_epsilon_two = (2._kp/cosh(x)**2)/mu**2

  end function sdi_epsilon_two


!epsilon_three(x)
  function sdi_epsilon_three(x,mu)    
    implicit none
    real(kp) :: sdi_epsilon_three
    real(kp), intent(in) :: x,mu
    
    sdi_epsilon_three = (-2._kp*tanh(x)**2)/mu**2
    
  end function sdi_epsilon_three

!returns x at the end of inflation
  function sdi_x_endinf(mu,xend)
    implicit none
    real(kp), intent(in) :: mu,xend
    real(kp) :: sdi_x_endinf
   
    sdi_x_endinf = xend
   
  end function sdi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function sdi_efold_primitive(x,mu)
    implicit none
    real(kp), intent(in) :: x,mu
    real(kp) :: sdi_efold_primitive

    if (mu.eq.0._kp) stop 'sdi_efold_primitive: mu=0!'

    sdi_efold_primitive = -mu**2*log(sinh(x))

  end function sdi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function sdi_x_trajectory(bfold,xend,mu)
    implicit none
    real(kp), intent(in) :: bfold, mu, xend
    real(kp) :: sdi_x_trajectory
    
    sdi_x_trajectory = asinh(sinh(xend)*exp(bfold/mu**2))
       
  end function sdi_x_trajectory


  
end module sdisr
