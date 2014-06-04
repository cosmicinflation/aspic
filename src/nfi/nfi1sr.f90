!slow-roll functions for the N-formalism inflation potential 1
!
!
!V(phi)) = M^4 exp (a x^b)
!
!1: a>0, b>1
!
!x = phi/Mp

module nfi1sr
  use infprec, only : kp
  use nficommon, only : nfi_norm_potential, nfi_norm_deriv_potential
  use nficommon, only : nfi_norm_deriv_second_potential
  use nficommon, only : nfi_epsilon_one, nfi_epsilon_two, nfi_epsilon_three
  use nficommon, only : nfi_x_epsoneunity, nfi_efold_primitive

  implicit none

  private

  public nfi1_norm_potential
  public nfi1_norm_deriv_potential, nfi1_norm_deriv_second_potential
  public nfi1_epsilon_one, nfi1_epsilon_two, nfi1_epsilon_three
  public nfi1_x_endinf, nfi1_efold_primitive

  public nfi1_numacc_xinimin

contains

  
!V/M^4
  function nfi1_norm_potential(x,a,b)
    implicit none
    real(kp) :: nfi1_norm_potential
    real(kp), intent(in) :: x,a,b

    nfi1_norm_potential = nfi_norm_potential(x,a,b)

  end function nfi1_norm_potential
 

!first derivative of the potential/M^4 with respect to x
  function nfi1_norm_deriv_potential(x,a,b)
    implicit none
    real(kp) :: nfi1_norm_deriv_potential
    real(kp), intent(in) :: x,a,b

    nfi1_norm_deriv_potential = nfi_norm_deriv_potential(x,a,b)

  end function nfi1_norm_deriv_potential


!second derivative of the potential/M^4 with respect to x
  function nfi1_norm_deriv_second_potential(x,a,b)
    implicit none
    real(kp) :: nfi1_norm_deriv_second_potential
    real(kp), intent(in) :: x,a,b

    nfi1_norm_deriv_second_potential = nfi_norm_deriv_second_potential(x,a,b)

  end function nfi1_norm_deriv_second_potential


!epsilon_one(x)
  function nfi1_epsilon_one(x,a,b)    
    implicit none
    real(kp) :: nfi1_epsilon_one
    real(kp), intent(in) :: x,a,b
    
    nfi1_epsilon_one = nfi_epsilon_one(x,a,b)
    
  end function nfi1_epsilon_one


!epsilon_two(x)
  function nfi1_epsilon_two(x,a,b)    
    implicit none
    real(kp) :: nfi1_epsilon_two
    real(kp), intent(in) :: x,a,b
    
    nfi1_epsilon_two = nfi_epsilon_two(x,a,b)
    
  end function nfi1_epsilon_two


!epsilon_three(x)
  function nfi1_epsilon_three(x,a,b)    
    implicit none
    real(kp) :: nfi1_epsilon_three
    real(kp), intent(in) :: x,a,b
    
    nfi1_epsilon_three = nfi_epsilon_three(x,a,b)
    
  end function nfi1_epsilon_three


!returns x and the end of inflation defined as epsilon1=1
  function nfi1_x_endinf(a,b)
    implicit none
    real(kp) :: nfi1_x_endinf
    real(kp), intent(in) :: a,b
    real(kp), dimension(2) :: xepsones

    xepsones = nfi_x_epsoneunity(a,b)

    nfi1_x_endinf = xepsones(1)

  end function nfi1_x_endinf


!return the minimal positive value of xini for ensuring eps1 >
!numerical accuracy
  function nfi1_numacc_xinimin(a,b)
    implicit none
    real(kp) :: nfi1_numacc_xinimin
    real(kp), intent(in) :: a,b

    if (b.eq.1._kp) stop 'nfi1_numacc_xinimin: b=1 is PLI'
    
    if (b.lt.1._kp) then
       nfi1_numacc_xinimin = 0._kp
    else
       nfi1_numacc_xinimin = (2._kp*epsilon(1._kp)/(a*a*b*b)) &
            ** (0.5_kp/(b-1._kp) )
    endif
    
  end function nfi1_numacc_xinimin



!this is integral[V(phi)/V'(phi) dphi]
  function nfi1_efold_primitive(x,a,b)
    implicit none
    real(kp), intent(in) :: x,a,b
    real(kp) :: nfi1_efold_primitive

    nfi1_efold_primitive = nfi_efold_primitive(x,a,b)

  end function nfi1_efold_primitive



!Return x as a function of bfold = N - N_end
  function nfi1_x_trajectory(bfold,a,b,xend)
    implicit none
    real(kp) :: nfi1_x_trajectory
    real(kp), intent(in) :: bfold,a,b,xend

    real(kp) :: efoldMax, xinimin

    xinimin = nfi1_numacc_xinimin(a,b)

    efoldMax = -nfi1_efold_primitive(xend,a,b) &
         +nfi1_efold_primitive(xinimin,a,b)

    if (-bfold.gt.efoldMax) then
       write(*,*)'nfi1_x_trajectory: not enough efolds!'
       write(*,*)'efold requested= efold maxi= ',-bfold,efoldMax
       if (b.gt.1._kp) write(*,*) 'xinimin (numacc)= ',xinimin
       stop 
    endif
    
    nfi1_x_trajectory = nfi_x_trajectory(bfold,a,b,xend)

  end function nfi1_x_trajectory


end module nfi1sr
