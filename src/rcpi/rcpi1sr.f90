!slow-roll function for radiatively corrected plateau inflation 1
!
!V(phi) = M^4 x^p [1 + alpha ln(x) + beta ln(x)^2]
!
!x=phi/Mp
!
!1: -V is positive definite (beta > 0 && alpha^2 < 4 beta)
!
!   -inflation is not transiently interrupted: x varies only in
!    well-defined slow-roll domains. It can be large only if inflation
!    is not retriggered at small field values, otherwise, it is
!    assumed to be confined in the small field domain. See function
!    rcpi1_xinimax(). Other cases cannot be studied within slow-roll.
!
module rcpi1sr
  use infprec, only : kp, tolkp
  use rcpicommon, only : rcpi_norm_potential, rcpi_norm_deriv_potential
  use rcpicommon, only : rcpi_norm_deriv_second_potential
  use rcpicommon, only : rcpi_epsilon_one, rcpi_epsilon_two, rcpi_epsilon_three
  use rcpicommon, only : rcpi_check_potzero, rcpi_check_derivpotzero
  use rcpicommon, only : rcpi_x_potzero, rcpi_x_derivpotzero, rcpi_efold_primitive
  use rcpicommon, only : rcpi_x_epsoneunity, rcpi_numacc_x_potbig

  implicit none

  private
  public rcpi1_check_params, rcpi1_efoldmax
  public rcpi1_norm_potential, rcpi1_norm_deriv_potential, rcpi1_norm_deriv_second_potential
  public rcpi1_epsilon_one, rcpi1_epsilon_two, rcpi1_epsilon_three
  public rcpi1_xinimax, rcpi1_x_endinf, rcpi1_x_trajectory  

contains


  function rcpi1_check_params(p,alpha,beta)
    implicit none
    logical :: rcpi1_check_params
    real(kp), intent(in) :: p, alpha, beta

    rcpi1_check_params = (beta.gt.0._kp) &
         .and.(alpha*alpha.le.4*beta) &
         .and.(p.ge.0._kp)

  end function rcpi1_check_params


  
  function rcpi1_norm_potential(x,p,alpha,beta)
    implicit none
    real(kp) :: rcpi1_norm_potential
    real(kp), intent(in) :: x,p,alpha,beta

    rcpi1_norm_potential = rcpi_norm_potential(x,p,alpha,beta)
    
  end function rcpi1_norm_potential

  

  function rcpi1_norm_deriv_potential(x,p,alpha,beta)
    implicit none
    real(kp) :: rcpi1_norm_deriv_potential
    real(kp), intent(in) :: x,p,alpha,beta

    rcpi1_norm_deriv_potential = rcpi_norm_deriv_potential(x,p,alpha,beta)
    
  end function rcpi1_norm_deriv_potential



  function rcpi1_norm_deriv_second_potential(x,p,alpha,beta)
    implicit none
    real(kp) :: rcpi1_norm_deriv_second_potential
    real(kp), intent(in) :: x,p,alpha,beta

    rcpi1_norm_deriv_second_potential = rcpi_norm_deriv_second_potential(x,p,alpha,beta)
    
  end function rcpi1_norm_deriv_second_potential

  

  function rcpi1_epsilon_one(x,p,alpha,beta)
    implicit none
    real(kp) :: rcpi1_epsilon_one
    real(kp), intent(in) :: x,p,alpha,beta

    rcpi1_epsilon_one = rcpi_epsilon_one(x,p,alpha,beta)

  end function rcpi1_epsilon_one


  
  function rcpi1_epsilon_two(x,p,alpha,beta)
    implicit none
    real(kp) :: rcpi1_epsilon_two
    real(kp), intent(in) :: x,p,alpha,beta

    rcpi1_epsilon_two = rcpi_epsilon_two(x,p,alpha,beta)

  end function rcpi1_epsilon_two


  
  function rcpi1_epsilon_three(x,p,alpha,beta)
    implicit none
    real(kp) :: rcpi1_epsilon_three
    real(kp), intent(in) :: x,p,alpha,beta

    rcpi1_epsilon_three = rcpi_epsilon_three(x,p,alpha,beta)

  end function rcpi1_epsilon_three
    

!returns the maximal value of x for inflation to not be interrupted
  function rcpi1_xinimax(p,alpha,beta)
    implicit none
    real(kp) :: rcpi1_xinimax
    real(kp), intent(in) :: p,alpha,beta
    real(kp), dimension(2) :: xdpotzero
    real(kp), dimension(5) :: xepsone

    real(kp), dimension(:), allocatable :: xeps

    integer :: neps, i,j

    if (.not.rcpi1_check_params(p,alpha,beta)) then
       stop 'rcpi1_xinimax: params out of range'
    endif

    if (rcpi_check_derivpotzero(p,alpha,beta)) then
!there is a false vacuum at xdpotzero(2) and false vacuum inflation
!could occurs by rolling down to it either from large field x>
!xdpotzero(2), or from the local maximum x > xdpotzero(1). We
!therefore remain confined to x < xdpotzero(1)
       xdpotzero = rcpi_x_derivpotzero(p,alpha,beta)
       rcpi1_xinimax = xdpotzero(1) - epsilon(1._kp)
       return
    endif

!the potential is monotonous increasing function of the field, but
!slow-roll could nevertheless still be interrupted by regions in which
!eps1>1

    xepsone = rcpi_x_epsoneunity(p,alpha,beta)
    neps = count(xepsone.ne.0._kp)
    allocate(xeps(neps))
    j=0
    do i=1,size(xepsone)
       if (xepsone(i).eq.0._kp) cycle
       j=j+1
       xeps(j) = xepsone(i)
    enddo

    select case(neps)

!no limit       
    case(1)
       rcpi1_xinimax = rcpi_numacc_x_potbig(p)

!large field inflation from x>>1 would be interrupted, small field
!inflation possible between xeps(1) and xeps(2)
    case(2,3)
       rcpi1_xinimax = xeps(2)
       
    case default
!for no potential local extremas, this should not happen
       stop 'rcpi1_xinimax: internal error!'
    end select

    deallocate(xeps)
        
  end function rcpi1_xinimax

  

  function rcpi1_x_endinf(p,alpha,beta)
    implicit none
    real(kp) :: rcpi1_x_endinf
    real(kp), intent(in) :: p, alpha, beta

!cf rcpi1_xinimax for more details
    real(kp), dimension(5) :: xepsone
    real(kp), dimension(:), allocatable :: xeps
    
    integer :: neps, i, j
    
    if (.not.rcpi1_check_params(p,alpha,beta)) then
       stop 'rcpi1_x_endinf: params out of range'
    endif

    xepsone = rcpi_x_epsoneunity(p,alpha,beta)
    neps = count(xepsone.ne.0._kp)
    allocate(xeps(neps))
    j=0
    do i=1,size(xepsone)
       if (xepsone(i).eq.0._kp) cycle
       j=j+1
       xeps(j) = xepsone(i)
    end do
       
    rcpi1_x_endinf = xeps(1)

    deallocate(xeps)
    
  end function rcpi1_x_endinf


  function rcpi1_efold_primitive(x,p,alpha,beta)
    implicit none
    real(kp) :: rcpi1_efold_primitive
    real(kp), intent(in) :: x,p,alpha,beta

    rcpi1_efold_primitive = rcpi_efold_primitive(x,p,alpha,beta)

  end function rcpi1_efold_primitive



  function rcpi1_x_trajectory(bfold,xend,p,alpha,beta)
    use infprec, only : transfert
    use inftools, only : zbrent
    use rcpicommon, only : find_rcpi_x_trajectory
    use rcpicommon, only : rcpi_efold_primitive
    implicit none
    real(kp) :: rcpi1_x_trajectory
    real(kp), intent(in) :: bfold,xend,p,alpha,beta
    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: rcpiData

    real(kp) :: mini, maxi

    mini = xend
    maxi = rcpi1_xinimax(p,alpha,beta)

    rcpiData%real1 = p
    rcpiData%real2 = alpha
    rcpiData%real3 = beta
    rcpiData%real4 = -bfold + rcpi_efold_primitive(xend,p,alpha,beta)

    rcpi1_x_trajectory = zbrent(find_rcpi_x_trajectory,mini,maxi,tolFind,rcpiData)

  end function rcpi1_x_trajectory

  

  function rcpi1_efoldmax(p,alpha,beta)
    implicit none
    real(kp) :: rcpi1_efoldmax
    real(kp), intent(in) :: p,alpha,beta

    real(kp) :: xinimax, xend

    xinimax = rcpi1_xinimax(p,alpha,beta)

    if (xinimax.eq.rcpi_numacc_x_potbig(p)) then
       rcpi1_efoldmax = huge(1._kp)
       return
    endif
    
    xend = rcpi1_x_endinf(p,alpha,beta)
    
    rcpi1_efoldmax = rcpi_efold_primitive(xinimax,p,alpha,beta) &
         - rcpi_efold_primitive(xend,p,alpha,beta)

  end function rcpi1_efoldmax
  
end module rcpi1sr
