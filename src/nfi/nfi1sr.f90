!slow-roll functions for the N-formalism inflation potential 1
!
!
!V(phi)) = M^4 exp (-a x^b)
!
!1: a>0, b>1
!
!x = phi/Mp

module nfi1sr
  use infprec, only : kp
  use nficommon, only : nfi_norm_potential, nfi_norm_deriv_potential
  use nficommon, only : nfi_norm_deriv_second_potential
  use nficommon, only : nfi_numacc_x_potbig, nfi_numacc_x_epsonenull
  use nficommon, only : nfi_epsilon_one, nfi_epsilon_two, nfi_epsilon_three
  use nficommon, only : nfi_x_epsoneunity, nfi_x_trajectory,nfi_efold_primitive
  use nficommon, only : logVbig
  
  implicit none

  private

  public nfi1_norm_potential
  public nfi1_norm_deriv_potential, nfi1_norm_deriv_second_potential
  public nfi1_epsilon_one, nfi1_epsilon_two, nfi1_epsilon_three
  public nfi1_x_endinf, nfi1_efold_primitive, nfi1_x_trajectory

  public nfi1_check_params, nfi1_numacc_amin, nfi1_amax
  public nfi1_numacc_xinimin, nfi1_numacc_amax

contains


  function nfi1_check_params(a,b)
    implicit none
    logical :: nfi1_check_params
    real(kp), intent(in) :: a,b

    nfi1_check_params = (a.gt.0._kp).and.(b.gt.1._kp)

  end function nfi1_check_params


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

    if (.not.nfi1_check_params(a,b)) then
       stop 'nfi1_x_endinf: nfi1 requires a>0, b>1'
    endif

    nfi1_x_endinf = nfi_x_epsoneunity(a,b)

  end function nfi1_x_endinf


!return the minimal positive value of xini for ensuring eps1 >
!numerical accuracy
  function nfi1_numacc_xinimin(a,b)
    implicit none
    real(kp) :: nfi1_numacc_xinimin
    real(kp), intent(in) :: a,b

    if (.not.nfi1_check_params(a,b)) then
       stop 'nfi1_numacc_xinimin: nfi1 requires a>0, b>1'
    endif

    nfi1_numacc_xinimin = nfi_numacc_x_epsonenull(a,b)


  end function nfi1_numacc_xinimin



!returns the minimal value of a given b such that xend < numacc_x_potbig
  function nfi1_numacc_amin(b)
    use nficommon, only : NfiBig
    implicit none
    real(kp) :: nfi1_numacc_amin
    real(kp), intent(in) :: b

    if (b.le.1._kp) then
       stop 'nfi1_numacc_amin: nfi1 requires b>1'
    endif

    nfi1_numacc_amin = log(NfiBig)**(1._kp-b)*(0.5_kp*b*b)**(-0.5_kp*b)
! Maximum 20 orders of magnitude in energy are allowed between
! 'star' and 'end', hence logVbig = 80 orders of magnitude for the
! potential
    nfi1_numacc_amin = logVbig*log(10._kp)**(1._kp-b)*(0.5_kp*b*b)**(-0.5_kp*b)

  end function nfi1_numacc_amin


!return the maximal value of a to get efold of inflation
!when b < 2 (for b>2, one can do an infinite number of efolds in x=0)
  function nfi1_amax(efold,b)
    implicit none
    real(kp) :: nfi1_amax
    real(kp), intent(in) :: efold,b
    logical, parameter :: display = .false.

    if (b.le.1._kp) then
       stop 'nfi1_amax: nfi1 requires b>1'
    endif

    if (b.ge.2._kp) then
       if (display) write(*,*)'nfi1_amax: for b>=2 amax is infinite!'
       nfi1_amax = huge(1._kp)
       return
    endif

    nfi1_amax = (b*(2._kp-b)*efold)**(1._kp-b) &
         * (0.5_kp*b*b)**(0.5_kp*(b-2._kp))


  end function nfi1_amax


!return the maximal value of a to get efold of inflation starting a
!xini=numacc (>0) to xend. Compared to nfi1_amax, this is due to
!numerical limitation
  function nfi1_numacc_amax(efold,b)
    use nficommon, only : NfiSmall
    implicit none
    real(kp) :: nfi1_numacc_amax
    real(kp), intent(in) :: efold,b

    if (b.le.1._kp) then
       stop 'nfi1_numacc_amax: nfi1 requires b>1'
    endif

    if (b.eq.2._kp) then
       nfi1_numacc_amax = -0.25*log(NfiSmall)/efold
    else
       nfi1_numacc_amax = (b*(2._kp-b)*efold &
            /(1._kp-NfiSmall**(0.5_kp*(2._kp-b)/(b-1._kp))) &
            )**(1._kp-b) &
            * (0.5_kp*b*b)**(0.5_kp*(b-2._kp))
    endif

  end function nfi1_numacc_amax


!this is integral[V(phi)/V'(phi) dphi]
  function nfi1_efold_primitive(x,a,b)
    implicit none
    real(kp), intent(in) :: x,a,b
    real(kp) :: nfi1_efold_primitive

    nfi1_efold_primitive = nfi_efold_primitive(x,a,b)

  end function nfi1_efold_primitive



!Return x as a function of bfold = N - N_end
  function nfi1_x_trajectory(bfold,xend,a,b)
    implicit none
    real(kp) :: nfi1_x_trajectory
    real(kp), intent(in) :: bfold,a,b,xend

    real(kp) :: efoldMax, xinimin

    if (.not.nfi1_check_params(a,b)) then
       stop 'nfi1_x_trajectory: nfi1 requires a>0, b>1'
    endif

    xinimin = nfi1_numacc_xinimin(a,b)

    efoldMax = -nfi1_efold_primitive(xend,a,b) &
         +nfi1_efold_primitive(xinimin,a,b)

    if (-bfold.gt.efoldMax) then
       write(*,*)'nfi1_x_trajectory: not enough efolds!'
       write(*,*)'efold requested= efold maxi= ',-bfold,efoldMax
       if (b.gt.1._kp) write(*,*) 'xinimin (numacc)= ',xinimin
       stop
    endif

    nfi1_x_trajectory = nfi_x_trajectory(bfold,xend,a,b)

  end function nfi1_x_trajectory


end module nfi1sr
