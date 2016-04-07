!equation of state inflation. Given 1+w(N)=1+P/rho, all formulae are
!exact and allow to get the field values and potential.

module eosflow
  use infprec, only : kp

  implicit none

  private

  public eos_x, eos_norm_potential
  public eos_norm_deriv_potential, eos_norm_deriv_second_potential
  public eos_epsilon_one, eos_epsilon_two, eos_epsilon_three
  

contains

!returns the negative branch of the of the field values (in reduced Planck units)
!given the primitive of sqrt[w(N)+1)] and up to a constant
  function eos_x(psqrwp1)
    implicit none
    real(kp) :: eos_x
    real(kp), intent(in) :: psqrwp1

    eos_x = -sqrt(3._kp)*psqrwp1

  end function eos_x


!returns the potential V(N) given the primitive of w(N)+1 and
!w(N)+1
  function eos_norm_potential(pwp1,wp1)
    implicit none
    real(kp) :: eos_norm_potential
    real(kp), intent(in) :: pwp1,wp1

    eos_norm_potential = (1._kp - 0.5_kp*wp1)*exp(-3._kp*pwp1)

  end function eos_norm_potential


!returns the absolute value of the derivative of the potential with respect to x (field
!values) up to a sign dV/dN * dx/dN
  function eos_norm_deriv_potential(pwp1,wp1,dwp1)
    implicit none
    real(kp) :: eos_norm_deriv_potential
    real(kp), intent(in) :: pwp1,wp1,dwp1

    real(kp) :: V

    V = eos_norm_potential(pwp1,wp1)

    eos_norm_deriv_potential = V/sqrt(3._kp*wp1) &
         * (3._kp*wp1 + 1._kp/(2._kp-wp1)*dwp1)

  end function eos_norm_deriv_potential


!returns the second derivative of the potential with respect to x
!(field values)
  function eos_norm_deriv_second_potential(pwp1,wp1,dwp1,d2wp1)
    implicit none
    real(kp) :: eos_norm_deriv_second_potential
    real(kp), intent(in) :: pwp1,wp1,dwp1,d2wp1
    
    real(kp) :: V

    V = eos_norm_potential(pwp1,wp1)

    eos_norm_deriv_second_potential = V/(3._kp*wp1) &
         * (9._kp*wp1*wp1 + 1.5_kp*(5._kp*wp1-2._kp)/(2._kp-wp1)*dwp1 &
         + 0.5_kp/(wp1*(2._kp-wp1))*dwp1*dwp1 - 1._kp/(2._kp-wp1)*d2wp1)

  end function eos_norm_deriv_second_potential



!first hubble flow function given w(N)+1
  function eos_epsilon_one(wp1)
    implicit none
    real(kp) :: eos_epsilon_one
    real(kp), intent(in) :: wp1

    eos_epsilon_one = 1.5_kp*wp1
    
  end function eos_epsilon_one


!second hubble flow function given w(N)+1 and d[w(N)+1]/dN
  function eos_epsilon_two(wp1,dwp1)
    implicit none
    real(kp) :: eos_epsilon_two
    real(kp), intent(in) :: wp1, dwp1

    eos_epsilon_two = dwp1/wp1

  end function eos_epsilon_two


!third hubble flow function given w(N)+1, d[w(N)+1]/dN and
!d^2[w(N)+1]/dN^2
  function eos_epsilon_three(wp1,dwp1,d2wp1)
    implicit none
    real(kp) :: eos_epsilon_three
    real(kp), intent(in) :: wp1,dwp1,d2wp1

    eos_epsilon_three = d2wp1/dwp1 - dwp1/wp1

  end function eos_epsilon_three



end module eosflow
