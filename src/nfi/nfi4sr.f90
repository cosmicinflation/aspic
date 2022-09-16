!slow-roll functions for the N-formalism inflation potential 4
!
!
!V(phi)) = M^4 exp (-a x^b)
!
!4: a>0, 0<b<1 and a<0, b<0
!
!x = phi/Mp

module nfi4sr
  use infprec, only : kp
  use nficommon, only : logVbig
  use nficommon, only : nfi_norm_potential, nfi_norm_deriv_potential
  use nficommon, only : nfi_norm_deriv_second_potential
  use nficommon, only : nfi_numacc_x_potbig, nfi_numacc_x_epsonenull
  use nficommon, only : nfi_epsilon_one, nfi_epsilon_two, nfi_epsilon_three
  use nficommon, only : nfi_x_epsoneunity, nfi_x_trajectory,nfi_efold_primitive

  implicit none

  private

  public nfi4_norm_potential
  public nfi4_norm_deriv_potential, nfi4_norm_deriv_second_potential
  public nfi4_epsilon_one, nfi4_epsilon_two, nfi4_epsilon_three
  public nfi4_efold_primitive, nfi4_x_trajectory

  public nfi4_check_params, nfi4_numacc_xendmax
  public nfi4_xinimin, nfi4_xendmin

contains


  function nfi4_check_params(a,b)
    implicit none
    logical :: nfi4_check_params
    real(kp), intent(in) :: a,b

    nfi4_check_params = ((a.gt.0._kp).and.(b.lt.1._kp).and.(b.gt.0._kp)) &
         .or. ((a.lt.0._kp).and.(b.lt.0._kp))

  end function nfi4_check_params


!V/M^4
  function nfi4_norm_potential(x,a,b)
    implicit none
    real(kp) :: nfi4_norm_potential
    real(kp), intent(in) :: x,a,b

    nfi4_norm_potential = nfi_norm_potential(x,a,b)

  end function nfi4_norm_potential


!first derivative of the potential/M^4 with respect to x
  function nfi4_norm_deriv_potential(x,a,b)
    implicit none
    real(kp) :: nfi4_norm_deriv_potential
    real(kp), intent(in) :: x,a,b

    nfi4_norm_deriv_potential = nfi_norm_deriv_potential(x,a,b)

  end function nfi4_norm_deriv_potential


!second derivative of the potential/M^4 with respect to x
  function nfi4_norm_deriv_second_potential(x,a,b)
    implicit none
    real(kp) :: nfi4_norm_deriv_second_potential
    real(kp), intent(in) :: x,a,b

    nfi4_norm_deriv_second_potential = nfi_norm_deriv_second_potential(x,a,b)

  end function nfi4_norm_deriv_second_potential


!epsilon_one(x)
  function nfi4_epsilon_one(x,a,b)
    implicit none
    real(kp) :: nfi4_epsilon_one
    real(kp), intent(in) :: x,a,b

    nfi4_epsilon_one = nfi_epsilon_one(x,a,b)

  end function nfi4_epsilon_one


!epsilon_two(x)
  function nfi4_epsilon_two(x,a,b)
    implicit none
    real(kp) :: nfi4_epsilon_two
    real(kp), intent(in) :: x,a,b

    nfi4_epsilon_two = nfi_epsilon_two(x,a,b)

  end function nfi4_epsilon_two


!epsilon_three(x)
  function nfi4_epsilon_three(x,a,b)
    implicit none
    real(kp) :: nfi4_epsilon_three
    real(kp), intent(in) :: x,a,b

    nfi4_epsilon_three = nfi_epsilon_three(x,a,b)

  end function nfi4_epsilon_three


!returns the minimal possible value of xini to start within the domain
!eps1<1
  function nfi4_xinimin(a,b)
    implicit none
    real(kp) :: nfi4_xinimin
    real(kp), intent(in) :: a,b

    nfi4_xinimin = nfi_x_epsoneunity(a,b)

  end function nfi4_xinimin



!returns the minimal possible value of xend to get efold number of
!inflation while still starting in the domain eps1<1, i.e. at xinimin
   function nfi4_xendmin(efold,a,b)
    implicit none
    real(kp) :: nfi4_xendmin
    real(kp), intent(in) :: efold,a,b

    real(kp) :: xinimin, xendmax, efoldMax

    if (.not.nfi4_check_params(a,b)) then
       stop 'nfi4_xendmax: nfi4 requires a>0, 0<b<1 or a<0, b<0'
    endif


    xinimin = nfi4_xinimin(a,b)
    xendmax = huge(1._kp)

    efoldMax = -nfi4_efold_primitive(xendmax,a,b) &
         + nfi4_efold_primitive(xinimin,a,b)

    if (efold.gt.efoldMax) then
       write(*,*)'nfi4_xendmin: not enough efolds!'
       write(*,*)'efold requested=',efold,'efold maxi= ',efoldMax
       stop
    endif

    nfi4_xendmin = nfi_x_trajectory(efold,xinimin,a,b)

  end function nfi4_xendmin


!returns the maximal value of xend such as eps(xend)>numacc and xend <
!xpotbig
  function nfi4_numacc_xendmax(a,b)
    implicit none
    real(kp) :: nfi4_numacc_xendmax
    real(kp), intent(in) :: a,b

    real(kp) :: xpotbig, xepsnull, xinimin

    if (.not.nfi4_check_params(a,b)) then
       stop 'nfi4_xendmax: nfi4 requires a>0, 0<b<1 or a<0, b<0'
    endif

    xpotbig = nfi_numacc_x_potbig(a,b)
    xepsnull = nfi_numacc_x_epsonenull(a,b)

!because if b<0, xepsnull->0 (the potential diverges in x=0)
    if (b.lt.0._kp) then
       nfi4_numacc_xendmax = xepsnull
    else
       nfi4_numacc_xendmax = min(xepsnull,xpotbig)
    endif

!check that there are no more than logVbig = 80 order of magnitude
!between Vmin and Vmax (since inflation must proceed between the
!Planck and BBN scales)
    xinimin = nfi4_xinimin(a,b)
    nfi4_numacc_xendmax = min (nfi4_numacc_xendmax, &
                      (xinimin**b+logVbig*log(10._kp)/a)**(1._kp/b))

  end function nfi4_numacc_xendmax



!this is integral[V(phi)/V'(phi) dphi]
  function nfi4_efold_primitive(x,a,b)
    implicit none
    real(kp), intent(in) :: x,a,b
    real(kp) :: nfi4_efold_primitive

    nfi4_efold_primitive = nfi_efold_primitive(x,a,b)

  end function nfi4_efold_primitive



!Return x as a function of bfold = N - N_end
  function nfi4_x_trajectory(bfold,xend,a,b)
    implicit none
    real(kp) :: nfi4_x_trajectory
    real(kp), intent(in) :: bfold,a,b,xend

    real(kp) :: efoldMax, xinimin

    if (.not.nfi4_check_params(a,b)) then
       stop 'nfi4_x_trajectory: nfi4 requires a>0, 0<b<1 or a<0, b<0'
    endif

    xinimin = nfi4_xinimin(a,b)

    efoldMax = -nfi4_efold_primitive(xend,a,b) &
         +nfi4_efold_primitive(xinimin,a,b)

    if (-bfold.gt.efoldMax) then
       write(*,*)'nfi4_x_trajectory: not enough efolds!'
       write(*,*)'efold requested= efold maxi= ',-bfold,efoldMax
       stop
    endif

    nfi4_x_trajectory = nfi_x_trajectory(bfold,xend,a,b)

  end function nfi4_x_trajectory


end module nfi4sr
