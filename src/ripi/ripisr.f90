!slow-roll functions for the renormalizable inflection point inflation potential
!
!V(phi) = M**4 (x**2 - x**3 + 9/32 x**4)
!
!x = phi/phi0

module ripisr
  use infprec, only : kp,pi,tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public  ripi_norm_potential, ripi_epsilon_one, ripi_epsilon_two, ripi_epsilon_three
  public  ripi_x_endinf, ripi_efold_primitive, ripi_x_trajectory, ripi_x_derivpotzero
  public  ripi_norm_deriv_potential, ripi_norm_deriv_second_potential

contains
  !returns V/M**4
  function ripi_norm_potential(x,phi0)
    implicit none
    real(kp) :: ripi_norm_potential
    real(kp), intent(in) :: x,phi0

    ripi_norm_potential = x**2-x**3+9._kp/32._kp*x**4

  end function ripi_norm_potential



  !returns the first derivative of the potential with respect to x, divided by M**4
  function ripi_norm_deriv_potential(x,phi0)
    implicit none
    real(kp) :: ripi_norm_deriv_potential
    real(kp), intent(in) :: x,phi0

    ripi_norm_deriv_potential = 1._kp/8._kp*x*(4._kp-3._kp*x)**2

  end function ripi_norm_deriv_potential



  !returns the *1._kp/cosond derivative of the potential with respect to x, divided by M**4
  function ripi_norm_deriv_second_potential(x,phi0)
    implicit none
    real(kp) :: ripi_norm_deriv_second_potential
    real(kp), intent(in) :: x,phi0

    ripi_norm_deriv_second_potential = 1._kp/8._kp*(-4._kp+3._kp*x)* &
         (-4._kp+9._kp*x)

  end function ripi_norm_deriv_second_potential



  !epsilon_one(x)
  function ripi_epsilon_one(x,phi0)    
    implicit none
    real(kp) :: ripi_epsilon_one
    real(kp), intent(in) :: x,phi0

    ripi_epsilon_one = ((8._kp*(4._kp-3._kp*x)**4)/(x**2*(32._kp+ &
         x*(-32._kp+9._kp*x))**2))/phi0**2

  end function ripi_epsilon_one


  !epsilon_two(x)
  function ripi_epsilon_two(x,phi0)    
    implicit none
    real(kp) :: ripi_epsilon_two
    real(kp), intent(in) :: x,phi0

    ripi_epsilon_two =((8._kp*(-4._kp+3._kp*x)*(-128._kp+x* &
         (160._kp+27._kp*x*(-4._kp+x))))/(x**2*(32._kp+ &
         x*(-32._kp+9._kp*x))**2))/phi0**2

  end function ripi_epsilon_two


  !epsilon_three(x)
  function ripi_epsilon_three(x,phi0)    
    implicit none
    real(kp) :: ripi_epsilon_three
    real(kp), intent(in) :: x,phi0

    ripi_epsilon_three = ((8._kp*(-4._kp+3._kp*x)*(16384._kp+3._kp* &
         x*(-16384._kp+x*(20992._kp+x* &
         (-15104._kp+27._kp*x*(256._kp+9._kp*x* &
         (-8._kp+x)))))))/(x**2*(32._kp+x* &
         (-32._kp+9._kp*x))**2*(-128._kp+x* &
         (160._kp+27._kp*x*(-4._kp+x)))))/phi0**2

  end function ripi_epsilon_three

!
  function ripi_x_derivpotzero(phi0)
    implicit none
    real(kp), intent(in), optional :: phi0
    real(kp) :: ripi_x_derivpotzero

    ripi_x_derivpotzero = 4._kp/3._kp

  end function ripi_x_derivpotzero

!inflection points
  function ripi_x_dderivpotzero(phi0)
    implicit none
    real(kp), intent(in), optional :: phi0
    real(kp), dimension(2) :: ripi_x_dderivpotzero

    ripi_x_dderivpotzero(1) = 4._kp/9._kp
    ripi_x_dderivpotzero(2) = 4._kp/3._kp

  end function ripi_x_dderivpotzero


  !returns x at the end of inflation defined as epsilon1=1
  function ripi_x_endinf(phi0)
    implicit none
    real(kp), intent(in) :: phi0
    real(kp) :: ripi_x_endinf

    ripi_x_endinf = real((1._kp/(27._kp*phi0))*2._kp*(9._kp*sqrt(2._kp) &
         +16._kp*phi0+(2._kp**(2._kp/3._kp)* &
         (-81._kp+18._kp*sqrt(2._kp)*phi0-20._kp* &
         phi0**2))/(phi0*(486._kp-297._kp* &
         sqrt(2._kp)*phi0+544._kp*phi0**2)+27._kp* &
         sqrt(2._kp)*(-27._kp+sqrt(3._kp)* &
         sqrt(cmplx(phi0**3*(-162._kp*sqrt(2._kp)+phi0* &
         (99._kp+64._kp*phi0*(-sqrt(2._kp)+phi0))),0._kp,kp))))** &
         (1._kp/3._kp)-2._kp**(1._kp/3._kp)*(phi0* &
         (486._kp-297._kp*sqrt(2._kp)*phi0+544._kp*phi0**2)+ &
         27._kp*sqrt(2._kp)*(-27._kp+sqrt(3._kp)* &
         sqrt(cmplx(phi0**3*(-162._kp*sqrt(2._kp)+phi0* &
         (99._kp+64._kp*phi0*(-sqrt(2._kp)+ phi0))),0._kp,kp)))) &
         **(1._kp/3._kp)))

  end function ripi_x_endinf



  !this is integral(V(phi)/V'(phi) dphi)
  function ripi_efold_primitive(x,phi0)
    implicit none
    real(kp), intent(in) :: x,phi0
    real(kp) :: ripi_efold_primitive

    if (x .eq. ripi_x_derivpotzero()) then
       stop 'ripi_efold_primitive evaluated at x=4/3 flat inflection point position!'
    endif

    ripi_efold_primitive = phi0**2*(-(2._kp*x/(9._kp)-x**2/8._kp+16._kp/(27._kp* &
         (3._kp*x-4._kp))+4._kp*log(4._kp-3._kp*x)/(27._kp)))

  end function ripi_efold_primitive


  !returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ripi_x_trajectory(bfold,xend,phi0)
    implicit none
    real(kp), intent(in) :: bfold, phi0, xend
    real(kp) :: ripi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ripiData


    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = ripi_x_derivpotzero()*(1._kp-epsilon(1._kp)) !Position of the flat inflection point

    ripiData%real1 = phi0
    ripiData%real2 = -bfold + ripi_efold_primitive(xend,phi0)

    ripi_x_trajectory = zbrent(find_ripi_x_trajectory,mini,maxi,tolFind,ripiData)

  end function ripi_x_trajectory

  function find_ripi_x_trajectory(x,ripiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ripiData
    real(kp) :: find_ripi_x_trajectory
    real(kp) :: phi0,NplusNuend

    phi0 = ripiData%real1
    NplusNuend = ripiData%real2

    find_ripi_x_trajectory = ripi_efold_primitive(x,phi0) - NplusNuend

  end function find_ripi_x_trajectory



end module ripisr
