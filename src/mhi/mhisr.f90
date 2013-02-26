!slow-roll functions for the mutated hilltop potential
!
!V(phi) = M^4 * [1 - sech(x)]
!
!x = phi/mu
!mu=mu/Mp

module mhisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  use inftools, only : zbrent
  implicit none

  private

  public  mhi_norm_potential, mhi_epsilon_one, mhi_epsilon_two, mhi_epsilon_three
  public  mhi_x_endinf, mhi_efold_primitive, mhi_x_trajectory
  public  mhi_norm_deriv_potential, mhi_norm_deriv_second_potential
 
contains
!returns V/M^4
  function mhi_norm_potential(x,mu)
    implicit none
    real(kp) :: mhi_norm_potential
    real(kp), intent(in) :: x
    real(kp), intent(in), optional :: mu

    mhi_norm_potential = 1._kp-1._kp/cosh(x)

  end function mhi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function mhi_norm_deriv_potential(x,mu)
    implicit none
    real(kp) :: mhi_norm_deriv_potential
    real(kp), intent(in) :: x
    real(kp), intent(in), optional :: mu

   mhi_norm_deriv_potential = tanh(x)/cosh(x)

  end function mhi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function mhi_norm_deriv_second_potential(x,mu)
    implicit none
    real(kp) :: mhi_norm_deriv_second_potential
    real(kp), intent(in) :: x
    real(kp), intent(in), optional :: mu

    mhi_norm_deriv_second_potential = (cosh(x)**(-3)-tanh(x)**2/cosh(x))

  end function mhi_norm_deriv_second_potential



!epsilon_one(x=phi/mu)
  function mhi_epsilon_one(x,mu)    
    implicit none
    real(kp) :: mhi_epsilon_one
    real(kp), intent(in) :: x,mu

  
    mhi_epsilon_one = 0.5_kp*(mu*tanh(x/2)*cosh(x))**(-2)

    
  end function mhi_epsilon_one


!epsilon_two(x=phi/mu)
  function mhi_epsilon_two(x,mu)    
    implicit none
    real(kp) :: mhi_epsilon_two
    real(kp), intent(in) :: x,mu
    
    mhi_epsilon_two = (sinh(x/2)**(-2)+2._kp*cosh(x)**(-2))/(mu**2)
    
  end function mhi_epsilon_two


!epsilon_three(x=phi/mu)
  function mhi_epsilon_three(x,mu)    
    implicit none
    real(kp) :: mhi_epsilon_three
    real(kp), intent(in) :: x,mu
    
    mhi_epsilon_three = (4._kp+cosh(x)+sinh(x/2)**(-2)-2._kp*cosh(x)**(-2))/ &
                        ((cosh(x)+sinh(x)**2)*mu**2)
    
  end function mhi_epsilon_three


!returns x=phi/mu at the end of inflation defined as epsilon1=1
  function mhi_x_endinf(mu)
    implicit none
    real(kp), intent(in) ::mu
    real(kp) :: mhi_x_endinf
    complex(kp) ::sechxEnd
    
    sechxEnd =  -1._kp/3._kp+(1._kp-6._kp*mu**2)/3._kp* &
                  (-1._kp+36._kp*mu**2+3._kp*sqrt(6._kp)*mu* &
                  sqrt(cmplx(4._kp*mu**4+22._kp*mu**2-1._kp,0._kp)))**(-1._kp/3._kp) &
                  +1._kp/3._kp*(-1._kp+36._kp*mu**2+3._kp*sqrt(6._kp)*mu* &
                  sqrt(cmplx(4._kp*mu**4+22._kp*mu**2-1._kp,0._kp)))**(1._kp/3._kp)

    mhi_x_endinf=acosh(1._kp/real(sechxEnd))

   
  end function mhi_x_endinf



!tmhis is integral[V(phi)/V'(phi) dphi]
  function mhi_efold_primitive(x,mu)
    implicit none
    real(kp), intent(in) :: x,mu
    real(kp) :: mhi_efold_primitive

    mhi_efold_primitive = mu**2*(cosh(x)-2._kp*log(cosh(x/2._kp)))

  end function mhi_efold_primitive


!returns x=phi/mi at bfold=-efolds before the end of inflation, ie N-Nend
  function mhi_x_trajectory(bfold,xend,mu)
    implicit none
    real(kp), intent(in) :: bfold, xend, mu
    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mhi_x_trajectory,mini,maxi
    type(transfert) :: mhiData

    mhiData%real1 = bfold
    mhiData%real2 = xend
    mhiData%real3 = mu

    mini = xend
    maxi = 100._kp
    
   ! mhi_x_trajectory = acosh(-1._kp-lambert(exp(-1._kp-cosh(xend))&
    !   *(1._kp+cosh(xend))*exp(bfold/(mu**2)),-1))  !Numerical Instabilities Using this analytical formula

    mhi_x_trajectory =zbrent(find_mhi_x_trajectory,mini,maxi,tolzbrent,mhiData)
       
  end function mhi_x_trajectory

  function find_mhi_x_trajectory(x,mhiData)   
    implicit none
    real(kp) :: find_mhi_x_trajectory
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: mhiData

    real(kp) :: bfold,xend,mu

    bfold = mhiData%real1
    xend = mhiData%real2
    mu = mhiData%real3

    find_mhi_x_trajectory = mhi_efold_primitive(x,mu) - mhi_efold_primitive(xend,mu) + bfold
  
  end function find_mhi_x_trajectory


end module mhisr
