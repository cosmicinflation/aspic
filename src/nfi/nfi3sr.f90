!slow-roll functions for the N-formalism inflation potential 3
!
!
!V(phi)) = M^4 exp (a x^b)
!
!3: a<0, 0<b<1 and a>0, b<0
!
!x = phi/Mp

module nfi3sr
  use infprec, only : kp
  use nficommon, only : nfi_norm_potential, nfi_norm_deriv_potential
  use nficommon, only : nfi_norm_deriv_second_potential
  use nficommon, only : nfi_numacc_x_potbig, nfi_numacc_x_epsonenull
  use nficommon, only : nfi_epsilon_one, nfi_epsilon_two, nfi_epsilon_three
  use nficommon, only : nfi_x_epsoneunity, nfi_x_trajectory,nfi_efold_primitive

  implicit none

  private

  public nfi3_norm_potential
  public nfi3_norm_deriv_potential, nfi3_norm_deriv_second_potential
  public nfi3_epsilon_one, nfi3_epsilon_two, nfi3_epsilon_three
  public nfi3_x_endinf, nfi3_efold_primitive, nfi3_x_trajectory

  public nfi3_check_params
  public nfi3_numacc_xinimax,  nfi3_numacc_absamax

contains


  function nfi3_check_params(a,b)
    implicit none
    logical :: nfi3_check_params
    real(kp), intent(in) :: a,b

    nfi3_check_params = ((a.lt.0._kp).and.(b.lt.1._kp).and.(b.gt.0._kp)) &
         .or. ((a.gt.0._kp).and.(b.lt.0._kp))

  end function nfi3_check_params

  
!V/M^4
  function nfi3_norm_potential(x,a,b)
    implicit none
    real(kp) :: nfi3_norm_potential
    real(kp), intent(in) :: x,a,b

    nfi3_norm_potential = nfi_norm_potential(x,a,b)

  end function nfi3_norm_potential
 

!first derivative of the potential/M^4 with respect to x
  function nfi3_norm_deriv_potential(x,a,b)
    implicit none
    real(kp) :: nfi3_norm_deriv_potential
    real(kp), intent(in) :: x,a,b

    nfi3_norm_deriv_potential = nfi_norm_deriv_potential(x,a,b)

  end function nfi3_norm_deriv_potential


!second derivative of the potential/M^4 with respect to x
  function nfi3_norm_deriv_second_potential(x,a,b)
    implicit none
    real(kp) :: nfi3_norm_deriv_second_potential
    real(kp), intent(in) :: x,a,b

    nfi3_norm_deriv_second_potential = nfi_norm_deriv_second_potential(x,a,b)

  end function nfi3_norm_deriv_second_potential


!epsilon_one(x)
  function nfi3_epsilon_one(x,a,b)    
    implicit none
    real(kp) :: nfi3_epsilon_one
    real(kp), intent(in) :: x,a,b
    
    nfi3_epsilon_one = nfi_epsilon_one(x,a,b)
    
  end function nfi3_epsilon_one


!epsilon_two(x)
  function nfi3_epsilon_two(x,a,b)    
    implicit none
    real(kp) :: nfi3_epsilon_two
    real(kp), intent(in) :: x,a,b
    
    nfi3_epsilon_two = nfi_epsilon_two(x,a,b)
    
  end function nfi3_epsilon_two


!epsilon_three(x)
  function nfi3_epsilon_three(x,a,b)    
    implicit none
    real(kp) :: nfi3_epsilon_three
    real(kp), intent(in) :: x,a,b
    
    nfi3_epsilon_three = nfi_epsilon_three(x,a,b)
    
  end function nfi3_epsilon_three


!returns x and the end of inflation defined as epsilon1=1
  function nfi3_x_endinf(a,b)
    implicit none
    real(kp) :: nfi3_x_endinf
    real(kp), intent(in) :: a,b

    if (.not.nfi3_check_params(a,b)) then
       stop 'nfi3_x_endinf: nfi3 requires a<0, 0<b<1 or a>0, b<0'
    endif

    nfi3_x_endinf = nfi_x_epsoneunity(a,b)

  end function nfi3_x_endinf


!return the minimal positive value of xini for ensuring eps1 >
!numerical accuracy
  function nfi3_numacc_xinimax(a,b)
    implicit none
    real(kp) :: nfi3_numacc_xinimax
    real(kp), intent(in) :: a,b

    if (.not.nfi3_check_params(a,b)) then
       stop 'nfi3_numacc_xinimax: nfi3 requires a<0, 0<b<1 or a>0, b<0'
    endif
      
    nfi3_numacc_xinimax = nfi_numacc_x_epsonenull(a,b)
       
  end function nfi3_numacc_xinimax



!returns the maximal value of |a| given b such that xinimax < numacc_x_potbig
  function nfi3_numacc_absamax(b)
    use nficommon, only : NfiBig
    implicit none
    real(kp) :: nfi3_numacc_absamax
    real(kp), intent(in) :: b

    if ((b.gt.1._kp).or.(b.eq.0._kp)) then
       stop 'nfi3_numacc_amin: nfi3 requires non-vanishing b<1'
    endif

    nfi3_numacc_absamax = log(NfiBig)**(1._kp-b) &
         *(0.5_kp*b*b/epsilon(1._kp))**(-0.5_kp*b)
    
  end function nfi3_numacc_absamax



!this is integral[V(phi)/V'(phi) dphi]
  function nfi3_efold_primitive(x,a,b)
    implicit none
    real(kp), intent(in) :: x,a,b
    real(kp) :: nfi3_efold_primitive

    nfi3_efold_primitive = nfi_efold_primitive(x,a,b)

  end function nfi3_efold_primitive



!Return x as a function of bfold = N - N_end
  function nfi3_x_trajectory(bfold,xend,a,b)
    implicit none
    real(kp) :: nfi3_x_trajectory
    real(kp), intent(in) :: bfold,a,b,xend

    real(kp) :: efoldMax, xinimax

    if (.not.nfi3_check_params(a,b)) then
       stop 'nfi3_x_trajectory: nfi3 requires a>0, b>1'
    endif

    xinimax = nfi3_numacc_xinimax(a,b)

    efoldMax = -nfi3_efold_primitive(xend,a,b) &
         +nfi3_efold_primitive(xinimax,a,b)

    if (-bfold.gt.efoldMax) then
       write(*,*)'nfi3_x_trajectory: not enough efolds!'
       write(*,*)'efold requested= efold maxi= ',-bfold,efoldMax
       stop 
    endif
    
    nfi3_x_trajectory = nfi_x_trajectory(bfold,xend,a,b)

  end function nfi3_x_trajectory


end module nfi3sr
