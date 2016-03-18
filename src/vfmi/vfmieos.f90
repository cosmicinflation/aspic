!equation of state functions for alpha-beta inflation. Combined with
!the module "common/eosflow.f90" and vfmireheat.f90, they can be used
!to compute all observable predictions in the Mukhanov's original
!parametric form.

module vfmieos

  use infprec, only : kp
  implicit none

  private

  public vfmi_eos, vfmi_deriv_eos, vfmi_deriv_second_eos
  public vfmi_primitive_eos, vfmi_primitive_sqrteos


contains

  

!returns w(bfold)+1 where bfold=-efold=N-Nend < 0
  function vfmi_eos(bfold,alpha,beta)
    implicit none
    real(kp) :: vfmi_eos
    real(kp), intent(in) :: bfold,alpha,beta

    vfmi_eos = beta/((1.5_kp*beta)**(1._kp/alpha) - bfold)**alpha
    
  end function vfmi_eos


!returns d(w+1)/dN
  function vfmi_deriv_eos(bfold,alpha,beta)
    implicit none
    real(kp) :: vfmi_deriv_eos
    real(kp), intent(in) :: bfold,alpha,beta

    vfmi_deriv_eos = alpha*beta &
         /((1.5_kp*beta)**(1._kp/alpha) - bfold)**(alpha+1._kp)

  end function vfmi_deriv_eos



!returns d^2(w+1)/dN^2
  function vfmi_deriv_second_eos(bfold,alpha,beta)
    implicit none
    real(kp) :: vfmi_deriv_second_eos
    real(kp), intent(in) :: bfold,alpha,beta

    vfmi_deriv_second_eos = alpha*(alpha+1._kp)*beta &
         /((1.5_kp*beta)**(1._kp/alpha) - bfold)**(alpha+2._kp)

  end function vfmi_deriv_second_eos


!returns the primitive \int (w+1) dN
  function vfmi_primitive_eos(bfold,alpha,beta)
    implicit none
    real(kp) :: vfmi_primitive_eos
    real(kp), intent(in) :: bfold,alpha,beta

    if (alpha.eq.1._kp) then
       vfmi_primitive_eos = -beta*log(1.5_kp*beta - bfold)
       return
    endif
       
    vfmi_primitive_eos = beta/(alpha-1._kp) &
         /((1.5_kp*beta)**(1._kp/alpha) - bfold)**(alpha-1._kp)

  end function vfmi_primitive_eos



!returns the primitive \int sqrt(w+1) dN
  function vfmi_primitive_sqrteos(bfold,alpha,beta)
    implicit none
    real(kp) :: vfmi_primitive_sqrteos
    real(kp), intent(in) :: bfold,alpha,beta

    if (alpha.eq.2._kp) then
       vfmi_primitive_sqrteos = -sqrt(beta)*log(sqrt(1.5_kp*beta) - bfold)
       return
    endif
       
    vfmi_primitive_sqrteos = -sqrt(beta)/(0.5_kp*alpha-1._kp) &
         /((1.5_kp*beta)**(1._kp/alpha) - bfold)**(0.5_kp*alpha-1._kp)
   

  end function vfmi_primitive_sqrteos


end module vfmieos
