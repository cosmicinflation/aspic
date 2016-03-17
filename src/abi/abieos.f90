!equation of state functions for alpha-beta inflation. Combined with
!the module "common/eosflow.f90" and abireheat.f90, they can be used
!to compute all observable predictions in the Mukhanov's original
!parametric form.

module abieos

  use infprec, only : kp
  implicit none

  private

  public abi_eos, abi_deriv_eos, abi_deriv_second_eos
  public abi_eos_primitive, abi_sqrteos_primitive


contains

  

!returns w(bfold)+1 where bfold=-efold=N-Nend < 0
  function abi_eos(bfold,alpha,beta)
    implicit none
    real(kp) :: abi_eos
    real(kp), intent(in) :: bfold,alpha,beta

    abi_eos = beta/((1.5_kp*beta)**(1._kp/alpha) - bfold)**alpha
    
  end function abi_eos


!returns d(w+1)/dN
  function abi_deriv_eos(bfold,alpha,beta)
    implicit none
    real(kp) :: abi_deriv_eos
    real(kp), intent(in) :: bfold,alpha,beta

    abi_deriv_eos = alpha*beta &
         /((1.5_kp*beta)**(1._kp/alpha) - bfold)**(alpha+1._kp)

  end function abi_deriv_eos



!returns d^2(w+1)/dN^2
  function abi_deriv_second_eos(bfold,alpha,beta)
    implicit none
    real(kp) :: abi_deriv_second_eos
    real(kp), intent(in) :: bfold,alpha,beta

    abi_deriv_second_eos = alpha*(alpha+1._kp)*beta &
         /((1.5_kp*beta)**(1._kp/alpha) - bfold)**(alpha+2._kp)

  end function abi_deriv_second_eos


!returns the primitive \int (w+1) dN
  function abi_eos_primitive(bfold,alpha,beta)
    implicit none
    real(kp) :: abi_eos_primitive
    real(kp), intent(in) :: bfold,alpha,beta

    if (alpha.eq.1._kp) then
       abi_eos_primitive = -beta*log(1.5_kp*beta - bfold)
       return
    endif
       
    abi_eos_primitive = beta/(alpha-1._kp) &
         /((1.5_kp*beta)**(1._kp/alpha) - bfold)**(alpha-1._kp)

  end function abi_eos_primitive



!returns the primitive \int sqrt(w+1) dN
  function abi_sqrteos_primitive(bfold,alpha,beta)
    implicit none
    real(kp) :: abi_sqrteos_primitive
    real(kp), intent(in) :: bfold,alpha,beta

    if (alpha.eq.2._kp) then
       abi_sqrteos_primitive = -sqrt(beta)*log(sqrt(1.5_kp*beta) - bfold)
       return
    endif
       
    abi_sqrteos_primitive = -sqrt(beta)/(0.5_kp*alpha-1._kp) &
         /((1.5_kp*beta)**(1._kp/alpha) - bfold)**(0.5_kp*alpha-1._kp)
   

  end function abi_sqrteos_primitive


end module abieos
