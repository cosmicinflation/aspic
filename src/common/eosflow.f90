module eosflow
  use infprec, only : kp

  implicit none

  private

  public eos_x, eos_norm_potential
  public eos_epsilon_one, eos_epsilon_two, eps_epsilon_three
  

contains

!returns the absolute value of the field (in reduced Planck units)
!given the primitive of sqrt[w(N)+1)] and up to a constant
  function eos_x(pSqrWp1)
    implicit none
    real(kp) :: eos_x
    real(kp), intent(in) :: pSqrWp1

    eos_x = sqrt(3._kp)*pSqrWp1

  end function eos_x


!returns the potential V(N) given the primitive of w(N)+1 and
!w(N)+1
  function eos_norm_potential(pWp1,wp1)
    implicit none
    real(kp) :: eos_norm_potential
    real(kp), intent(in) :: pwp1,wp1

    eos_norm_potential = (1._kp - 0.5_kp*wp1)*exp(-3._kp*pWp1)

  end function eos_norm_potential


!first hubble flow function given w(N)+1
  function eos_epsilon_one(wp1)
    implicit none
    real(kp) :: eos_epsilon_one
    real(kp), intent(in) :: wp1

    eos_epsilon_one = 1.5_kp*wp1
    
  end function eos_epsilon_one


!second hubble flow function given w(N)+1 and d[w(N)+1]/dN
  function eos_epsilon_two(wp1,dWp1)
    implicit none
    real(kp) :: eos_epsilon_two
    real(kp), intent(in) :: wp1, dWp1

    eos_epsilon_two = dWp1/wp1

  end function eos_epsilon_two


!third hubble flow function given w(N)+1, d[w(N)+1]/dN and
!d^2[w(N)+1]/dN^2
  function eos_epsilon_three(wp1,dWp1,d2Wp1)
    implicit none
    real(kp) :: eos_epsilon_three
    real(kp), intent(in) :: wp1,dWp1,d2Wp1

    eos_epsilon_three = d2Wp1/dWp1 - dWp1/wp1

  end function eos_epsilon_three



end module eosflow
