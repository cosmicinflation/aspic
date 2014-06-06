!slow-roll functions for the N-formalism inflation potential 2
!
!
!V(phi)) = M^4 exp (a x^b)
!
!2: a<0, b>1
!
!x = phi/Mp

module nfi2sr
  use infprec, only : kp
  use nficommon, only : nfi_norm_potential, nfi_norm_deriv_potential
  use nficommon, only : nfi_norm_deriv_second_potential
  use nficommon, only : nfi_numacc_x_potbig, nfi_numacc_x_epsonenull
  use nficommon, only : nfi_epsilon_one, nfi_epsilon_two, nfi_epsilon_three
  use nficommon, only : nfi_x_epsoneunity, nfi_x_trajectory,nfi_efold_primitive

  implicit none

  private

  public nfi2_norm_potential
  public nfi2_norm_deriv_potential, nfi2_norm_deriv_second_potential
  public nfi2_epsilon_one, nfi2_epsilon_two, nfi2_epsilon_three
  public nfi2_efold_primitive, nfi2_x_trajectory

  public nfi2_check_params
  public nfi2_numacc_xinimax, nfi2_numacc_xendmax, nfi2_numacc_xendmin
  public nfi2_xinimax, nfi2_xendmax, nfi2_amin

contains


  function nfi2_check_params(a,b)
    implicit none
    logical :: nfi2_check_params
    real(kp), intent(in) :: a,b

    nfi2_check_params = (a.lt.0._kp).and.(b.gt.1._kp)

  end function nfi2_check_params

  
!V/M^4
  function nfi2_norm_potential(x,a,b)
    implicit none
    real(kp) :: nfi2_norm_potential
    real(kp), intent(in) :: x,a,b

    nfi2_norm_potential = nfi_norm_potential(x,a,b)

  end function nfi2_norm_potential
 

!first derivative of the potential/M^4 with respect to x
  function nfi2_norm_deriv_potential(x,a,b)
    implicit none
    real(kp) :: nfi2_norm_deriv_potential
    real(kp), intent(in) :: x,a,b

    nfi2_norm_deriv_potential = nfi_norm_deriv_potential(x,a,b)

  end function nfi2_norm_deriv_potential


!second derivative of the potential/M^4 with respect to x
  function nfi2_norm_deriv_second_potential(x,a,b)
    implicit none
    real(kp) :: nfi2_norm_deriv_second_potential
    real(kp), intent(in) :: x,a,b

    nfi2_norm_deriv_second_potential = nfi_norm_deriv_second_potential(x,a,b)

  end function nfi2_norm_deriv_second_potential


!epsilon_one(x)
  function nfi2_epsilon_one(x,a,b)    
    implicit none
    real(kp) :: nfi2_epsilon_one
    real(kp), intent(in) :: x,a,b
    
    nfi2_epsilon_one = nfi_epsilon_one(x,a,b)
    
  end function nfi2_epsilon_one


!epsilon_two(x)
  function nfi2_epsilon_two(x,a,b)    
    implicit none
    real(kp) :: nfi2_epsilon_two
    real(kp), intent(in) :: x,a,b
    
    nfi2_epsilon_two = nfi_epsilon_two(x,a,b)
    
  end function nfi2_epsilon_two


!epsilon_three(x)
  function nfi2_epsilon_three(x,a,b)    
    implicit none
    real(kp) :: nfi2_epsilon_three
    real(kp), intent(in) :: x,a,b
    
    nfi2_epsilon_three = nfi_epsilon_three(x,a,b)
    
  end function nfi2_epsilon_three


!returns the maximal possible value of xini to start within the domain
!eps1<1
  function nfi2_xinimax(a,b)
    implicit none
    real(kp) :: nfi2_xinimax
    real(kp), intent(in) :: a,b

    nfi2_xinimax = nfi_x_epsoneunity(a,b)
    
  end function nfi2_xinimax



!returns the maximal possible value of xend to get efold number of
!inflation while still starting in the domain eps1<1
   function nfi2_xendmax(efold,a,b)
    implicit none
    real(kp) :: nfi2_xendmax
    real(kp), intent(in) :: efold,a,b
    
    real(kp) :: xinimax, xendmin, efoldMax

    if (.not.nfi2_check_params(a,b)) then
       stop 'nfi2_xendmax: nfi2 requires a<0, b>1'
    endif


    xinimax = nfi2_xinimax(a,b)
    xendmin = tiny(0._kp)

    efoldMax = -nfi2_efold_primitive(xendmin,a,b) &
         + nfi2_efold_primitive(xinimax,a,b)
    
    if (efold.gt.efoldMax) then
       write(*,*)'nfi2_xendmax: not enough efolds!'
       write(*,*)'efold requested=',efold,'efold maxi= ',efoldMax
       stop
    endif

    nfi2_xendmax = nfi_x_trajectory(efold,xinimax,a,b)

  end function nfi2_xendmax



!return the maximal positive value of xini both inflation and numerical
!value for the potential (< huge)
  function nfi2_numacc_xinimax(a,b)
    implicit none
    real(kp) :: nfi2_numacc_xinimax
    real(kp), intent(in) :: a,b

    if (.not.nfi2_check_params(a,b)) then
       stop 'nfi1_numacc_xinimax: nfi2 requires a<0, b>1'
    endif

    nfi2_numacc_xinimax = min(nfi2_xinimax(a,b),nfi_numacc_x_potbig(a,b))
   
  end function nfi2_numacc_xinimax



!returns the maximal possible value of xend to get efold number of
!inflation while still starting in the domain xini < numacc_xinimax
   function nfi2_numacc_xendmax(efold,a,b)
    implicit none
    real(kp) :: nfi2_numacc_xendmax
    real(kp), intent(in) :: efold,a,b
    
    real(kp) :: xinimax, xendmin, efoldMax

    if (.not.nfi2_check_params(a,b)) then
       stop 'nfi2_numacc_xendmax: nfi2 requires a<0, b>1'
    endif


    xinimax = nfi2_numacc_xinimax(a,b)
    xendmin = tiny(0._kp)

    efoldMax = -nfi2_efold_primitive(xendmin,a,b) &
         + nfi2_efold_primitive(xinimax,a,b)
    
    if (efold.gt.efoldMax) then
       write(*,*)'nfi2_numacc_xendmax: not enough efolds!'
       write(*,*)'efold requested=',efold,'efold maxi= ',efoldMax
       stop
    endif

    nfi2_numacc_xendmax = nfi_x_trajectory(efold,xinimax,a,b)

  end function nfi2_numacc_xendmax




!return the minimal positive value of xend for ensuring numerical
!value of eps1>numprec
  function nfi2_numacc_xendmin(a,b)
    implicit none
    real(kp) :: nfi2_numacc_xendmin
    real(kp), intent(in) :: a,b

    if (.not.nfi2_check_params(a,b)) then
       stop 'nfi1_numacc_xendmin: nfi2 requires a<0, b>1'
    endif

    nfi2_numacc_xendmin = nfi_numacc_x_epsonenull(a,b)
   
  end function nfi2_numacc_xendmin



!return the minimal value of a to get efoldNum efolds of inflation
!when b < 2 (for b>2, one can do an infinite number of efolds by
!pushing xend close to 0)
  function nfi2_amin(efold,b)
    implicit none
    real(kp) :: nfi2_amin
    real(kp), intent(in) :: efold,b

    if (b.le.1._kp) then
       stop 'nfi2_amin: nfi2 requires b>1'
    endif
    
    if (b.ge.2._kp) then
       write(*,*)'nfi2_amin: for b>=2 amin is -infinite!'
       nfi2_amin = -huge(1._kp)
       return
    endif

    nfi2_amin = -(b*(2._kp-b)*efold)**(1._kp-b) &
         * (0.5_kp*b*b)**(0.5_kp*(b-2._kp))


  end function nfi2_amin


!this is integral[V(phi)/V'(phi) dphi]
  function nfi2_efold_primitive(x,a,b)
    implicit none
    real(kp), intent(in) :: x,a,b
    real(kp) :: nfi2_efold_primitive

    nfi2_efold_primitive = nfi_efold_primitive(x,a,b)

  end function nfi2_efold_primitive



!Return x as a function of bfold = N - N_end
  function nfi2_x_trajectory(bfold,xend,a,b)
    implicit none
    real(kp) :: nfi2_x_trajectory
    real(kp), intent(in) :: bfold,a,b,xend

    real(kp) :: efoldMax, xinimax

    if (.not.nfi2_check_params(a,b)) then
       stop 'nfi2_x_trajectory: nfi2 requires a<0, b>1'
    endif
    

    xinimax = nfi2_xinimax(a,b)

    efoldMax = -nfi2_efold_primitive(xend,a,b) &
         +nfi2_efold_primitive(xinimax,a,b)

    if (-bfold.gt.efoldMax) then
       write(*,*)'nfi2_x_trajectory: not enough efolds!'
       write(*,*)'efold requested= efold maxi= ',-bfold,efoldMax
       stop
    endif
    
    nfi2_x_trajectory = nfi_x_trajectory(bfold,xend,a,b)

  end function nfi2_x_trajectory


end module nfi2sr
