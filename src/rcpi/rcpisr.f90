!slow-roll function for radiatively corrected plateau inflation
!
!V(phi) = M^4 x^p [1 + alpha ln(x) + beta ln(x)^2]
!
!x=phi/Mp
!
!V is positive definite (beta > 0 && alpha^2 < 4 beta)
!
!inflation is not transiently interrupted: x varies only in
!well-defined slow-roll domains. It can be large only if inflation
!is not retriggered at small field values, otherwise, it is
!assumed to be confined in the small field domain. See function
!rcpi_xinimax(). Other cases cannot be studied within slow-roll.
!
module rcpisr
  use infprec, only : kp, tolkp
  use rcpicommon, only : rcpi_norm_potential, rcpi_norm_deriv_potential
  use rcpicommon, only : rcpi_norm_deriv_second_potential
  use rcpicommon, only : rcpi_epsilon_one, rcpi_epsilon_two, rcpi_epsilon_three
  use rcpicommon, only : rcpi_check_potzero, rcpi_check_derivpotzero
  use rcpicommon, only : rcpi_x_potzero, rcpi_x_derivpotzero, rcpi_efold_primitive
  use rcpicommon, only : rcpi_x_epsoneunity, rcpi_numacc_x_potbig

  
  implicit none

  real(kp), parameter :: potzero = 100._kp*epsilon(1._kp)

  private
  public rcpi_check_params, rcpi_efoldmax
  public rcpi_norm_potential, rcpi_norm_deriv_potential, rcpi_norm_deriv_second_potential
  public rcpi_epsilon_one, rcpi_epsilon_two, rcpi_epsilon_three
  public rcpi_xinimax, rcpi_x_endinf, rcpi_x_trajectory, rcpi_efold_primitive
  public rcpi_numacc_alphamax, rcpi_numacc_alphamin, rcpi_numacc_betamin

contains


  function rcpi_check_params(p,alpha,beta)
    implicit none
    logical :: rcpi_check_params
    real(kp), intent(in) :: p, alpha, beta

    rcpi_check_params = (beta.gt.0._kp) &
         .and.(abs(alpha).le.2._kp*sqrt(beta)) &
         .and.(p.ge.0._kp)

  end function rcpi_check_params


!return the maximal value of alpha < \sqrt(beta) such that the minimal
!value of the potential (very close to zero) can be numerically
!tracked down. Above this limit, quantities such as V'/V are no longer
!properly calculated...
  function rcpi_numacc_alphamax(p,beta)
    implicit none
    real(kp) :: rcpi_numacc_alphamax
    real(kp), intent(in) :: p,beta
    real(kp) :: epsnumacc

    if (beta.eq.0._kp) stop 'rcpi_numacc_alphamax: beta=0!'
    
    epsnumacc = min(max(1._kp,exp(p/sqrt(beta))) * potzero,0.5_kp)
    
    rcpi_numacc_alphamax = 2._kp*sqrt(beta*(1._kp - epsnumacc))

  end function rcpi_numacc_alphamax

!same as above for minimal values (non-symmetric!!)
  function rcpi_numacc_alphamin(p,beta)
    implicit none
    real(kp) :: rcpi_numacc_alphamin
    real(kp), intent(in) :: p,beta
    real(kp) :: epsnumacc
    
    if (beta.eq.0._kp) stop 'rcpi_numacc_alphamin: beta=0!'

    epsnumacc = max(1._kp,exp(-p/sqrt(beta))) * potzero
    
    rcpi_numacc_alphamin = -2._kp*sqrt(beta*(1._kp - epsnumacc))

  end function rcpi_numacc_alphamin
  
!for alpha of order 1, the potential becomes very close to zero at its
!extremal value for small values of beta
  function rcpi_numacc_betamin(p)
    implicit none
    real(kp) :: rcpi_numacc_betamin
    real(kp), intent(in) :: p

    rcpi_numacc_betamin = -p/log(potzero)

  end function rcpi_numacc_betamin

  
!returns the maximal value of x for inflation to not be interrupted
  function rcpi_xinimax(p,alpha,beta)
    implicit none
    real(kp) :: rcpi_xinimax
    real(kp), intent(in) :: p,alpha,beta
    real(kp), dimension(2) :: xdpotzero
    real(kp), dimension(5) :: xepsone

    real(kp), dimension(:), allocatable :: xeps

    integer :: neps, i,j

    if (.not.rcpi_check_params(p,alpha,beta)) then
       stop 'rcpi_xinimax: params out of range'
    endif

    if (rcpi_check_derivpotzero(p,alpha,beta)) then
!there is a false vacuum at xdpotzero(2) and false vacuum inflation
!could occurs by rolling down to it either from large field x>
!xdpotzero(2), or from the local maximum x > xdpotzero(1). We
!therefore remain confined to x < xdpotzero(1)
       xdpotzero = rcpi_x_derivpotzero(p,alpha,beta)
       rcpi_xinimax = xdpotzero(1)*(1._kp - epsilon(1._kp))
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
       rcpi_xinimax = rcpi_numacc_x_potbig(p)

!large field inflation from x>>1 would be interrupted, small field
!inflation possible between xeps(1) and xeps(2)
    case(2,3)
       rcpi_xinimax = xeps(2)
       
    case default
!for no potential local extremas, this should not happen
       stop 'rcpi_xinimax: internal error!'
    end select
    
    deallocate(xeps)
        
  end function rcpi_xinimax

  

  function rcpi_x_endinf(p,alpha,beta)
    implicit none
    real(kp) :: rcpi_x_endinf
    real(kp), intent(in) :: p, alpha, beta

!cf rcpi_xinimax for more details
    real(kp), dimension(5) :: xepsone
    real(kp), dimension(:), allocatable :: xeps
    
    integer :: neps, i, j
    
    if (.not.rcpi_check_params(p,alpha,beta)) then
       write(*,*)'p= alpha= beta= ',p,alpha,beta
       stop 'rcpi_x_endinf: params out of range'
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
       
    rcpi_x_endinf = xeps(1)

    deallocate(xeps)
    
  end function rcpi_x_endinf




  function rcpi_x_trajectory(bfold,xend,p,alpha,beta)
    use infprec, only : transfert
    use inftools, only : zbrent
    use rcpicommon, only : find_rcpi_x_trajectory
    use rcpicommon, only : rcpi_efold_primitive
    implicit none
    real(kp) :: rcpi_x_trajectory
    real(kp), intent(in) :: bfold,xend,p,alpha,beta
    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: rcpiData

    real(kp) :: mini, maxi

    mini = xend
    maxi = rcpi_xinimax(p,alpha,beta)

    rcpiData%real1 = p
    rcpiData%real2 = alpha
    rcpiData%real3 = beta
    rcpiData%real4 = -bfold + rcpi_efold_primitive(xend,p,alpha,beta)

    rcpi_x_trajectory = zbrent(find_rcpi_x_trajectory,mini,maxi,tolFind,rcpiData)

  end function rcpi_x_trajectory

  

  function rcpi_efoldmax(p,alpha,beta)
    implicit none
    real(kp) :: rcpi_efoldmax
    real(kp), intent(in) :: p,alpha,beta

    real(kp) :: xinimax, xend

    xinimax = rcpi_xinimax(p,alpha,beta)

    if (xinimax.eq.rcpi_numacc_x_potbig(p)) then
       rcpi_efoldmax = huge(1._kp)
       return
    endif
    
    xend = rcpi_x_endinf(p,alpha,beta)

    rcpi_efoldmax = rcpi_efold_primitive(xinimax,p,alpha,beta) &
         - rcpi_efold_primitive(xend,p,alpha,beta)

  end function rcpi_efoldmax
  
end module rcpisr
