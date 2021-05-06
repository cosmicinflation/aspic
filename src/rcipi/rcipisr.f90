!slow-roll function for radiatively corrected inflection point inflation
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
!rcipi_xinimax(). Other cases cannot be studied within slow-roll.
!
module rcipisr
  use infprec, only : kp, tolkp
  use rcipicommon, only : rcipi_norm_potential, rcipi_norm_deriv_potential
  use rcipicommon, only : rcipi_norm_deriv_second_potential
  use rcipicommon, only : rcipi_epsilon_one, rcipi_epsilon_two, rcipi_epsilon_three
  use rcipicommon, only : rcipi_check_potzero, rcipi_check_derivpotzero
  use rcipicommon, only : rcipi_x_potzero, rcipi_x_derivpotzero, rcipi_efold_primitive
  use rcipicommon, only : rcipi_x_epsoneunity, rcipi_numacc_x_potbig
  use rcipicommon, only : rcipi_alpha_zero

  
  implicit none

  real(kp), parameter :: potzero = 100._kp*epsilon(1._kp)

  private
  public rcipi_check_params, rcipi_efoldmax
  public rcipi_norm_potential, rcipi_norm_deriv_potential, rcipi_norm_deriv_second_potential
  public rcipi_epsilon_one, rcipi_epsilon_two, rcipi_epsilon_three, rcipi_alpha_zero
  public rcipi_xinimax, rcipi_x_endinf, rcipi_x_trajectory, rcipi_efold_primitive
  public rcipi_numacc_alphamax, rcipi_numacc_alphamin, rcipi_numacc_betamin

contains


  function rcipi_check_params(p,alpha,beta)
    implicit none
    logical :: rcipi_check_params
    real(kp), intent(in) :: p, alpha, beta

    rcipi_check_params = (beta.gt.0._kp) &
         .and.(abs(alpha).le.2._kp*sqrt(beta)) &
         .and.(p.ge.0._kp)

  end function rcipi_check_params


!return the maximal value of alpha < \sqrt(beta) such that the minimal
!value of the potential (very close to zero) can be numerically
!tracked down. Above this limit, quantities such as V'/V are no longer
!properly calculated...
  function rcipi_numacc_alphamax(p,beta)
    implicit none
    real(kp) :: rcipi_numacc_alphamax
    real(kp), intent(in) :: p,beta
    real(kp) :: epsnumacc

    if (beta.eq.0._kp) stop 'rcipi_numacc_alphamax: beta=0!'
    
    epsnumacc = min(max(1._kp,exp(p/sqrt(beta))) * potzero,0.5_kp)
    
    rcipi_numacc_alphamax = 2._kp*sqrt(beta*(1._kp - epsnumacc))

  end function rcipi_numacc_alphamax

!same as above for minimal values (non-symmetric!!)
  function rcipi_numacc_alphamin(p,beta)
    implicit none
    real(kp) :: rcipi_numacc_alphamin
    real(kp), intent(in) :: p,beta
    real(kp) :: epsnumacc
    
    if (beta.eq.0._kp) stop 'rcipi_numacc_alphamin: beta=0!'

    epsnumacc = max(1._kp,exp(-p/sqrt(beta))) * potzero
    
    rcipi_numacc_alphamin = -2._kp*sqrt(beta*(1._kp - epsnumacc))

  end function rcipi_numacc_alphamin
  
!for alpha of order 1, the potential becomes very close to zero at its
!extremal value for small values of beta
  function rcipi_numacc_betamin(p)
    implicit none
    real(kp) :: rcipi_numacc_betamin
    real(kp), intent(in) :: p

    rcipi_numacc_betamin = -p/log(potzero)

  end function rcipi_numacc_betamin

  
!returns the maximal value of x for inflation to not be interrupted
  function rcipi_xinimax(p,alpha,beta)
    implicit none
    real(kp) :: rcipi_xinimax
    real(kp), intent(in) :: p,alpha,beta
    real(kp), dimension(2) :: xdpotzero
    real(kp), dimension(5) :: xepsone

    real(kp), dimension(:), allocatable :: xeps

    integer :: neps, i,j

    if (.not.rcipi_check_params(p,alpha,beta)) then
       write(*,*)'p= alpha= beta= ',p,alpha,beta
       stop 'rcipi_xinimax: params out of range'
    endif

    if (rcipi_check_derivpotzero(p,alpha,beta)) then
!there is a false vacuum at xdpotzero(2) and false vacuum inflation
!could occurs by rolling down to it either from large field x>
!xdpotzero(2), or from the local maximum x > xdpotzero(1). We
!therefore remain confined to x < xdpotzero(1)
       xdpotzero = rcipi_x_derivpotzero(p,alpha,beta)
       rcipi_xinimax = xdpotzero(1)*(1._kp - epsilon(1._kp))
!       rcipi_xinimax = xdpotzero(1)- epsilon(1._kp)
       return
    endif

!the potential is monotonous increasing function of the field, but
!slow-roll could nevertheless still be interrupted by regions in which
!eps1>1

    xepsone = rcipi_x_epsoneunity(p,alpha,beta)
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
       rcipi_xinimax = rcipi_numacc_x_potbig(p)

!large field inflation from x>>1 is interrupted, but towards the end
    case(2,3)
       rcipi_xinimax = rcipi_numacc_x_potbig(p)
       
    case default
!for no potential local extremas, this should not happen
       stop 'rcipi_xinimax: internal error!'
    end select
    
    deallocate(xeps)
        
  end function rcipi_xinimax

  

  function rcipi_x_endinf(p,alpha,beta)
    implicit none
    real(kp) :: rcipi_x_endinf
    real(kp), intent(in) :: p, alpha, beta

!cf rcipi_xinimax for more details
    real(kp), dimension(5) :: xepsone
    real(kp), dimension(:), allocatable :: xeps
    
    integer :: neps, i, j
    
    if (.not.rcipi_check_params(p,alpha,beta)) then
       write(*,*)'p= alpha= beta= ',p,alpha,beta
       stop 'rcipi_x_endinf: params out of range'
    endif

    xepsone = rcipi_x_epsoneunity(p,alpha,beta)
    
    neps = count(xepsone.ne.0._kp)
    allocate(xeps(neps))
    j=0
    do i=1,size(xepsone)
       if (xepsone(i).eq.0._kp) cycle
       j=j+1
       xeps(j) = xepsone(i)
    end do

!we return the smallest root (ordered done in rcipi_x_epsoneunity)
    rcipi_x_endinf = xeps(1)

    deallocate(xeps)
    
  end function rcipi_x_endinf




  function rcipi_x_trajectory(bfold,xend,p,alpha,beta)
    use infprec, only : transfert
    use inftools, only : zbrent
    use rcipicommon, only : find_rcipi_x_trajectory
    use rcipicommon, only : rcipi_efold_primitive
    implicit none
    real(kp) :: rcipi_x_trajectory
    real(kp), intent(in) :: bfold,xend,p,alpha,beta
    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: rcipiData

    real(kp) :: mini, maxi

    mini = xend
    maxi = rcipi_xinimax(p,alpha,beta)

    rcipiData%real1 = p
    rcipiData%real2 = alpha
    rcipiData%real3 = beta
    rcipiData%real4 = -bfold + rcipi_efold_primitive(xend,p,alpha,beta)

    rcipi_x_trajectory = zbrent(find_rcipi_x_trajectory,mini,maxi,tolFind,rcipiData)

  end function rcipi_x_trajectory

  

  function rcipi_efoldmax(p,alpha,beta)
    implicit none
    real(kp) :: rcipi_efoldmax
    real(kp), intent(in) :: p,alpha,beta

    real(kp) :: xinimax, xend

    xinimax = rcipi_xinimax(p,alpha,beta)

    if (xinimax.eq.rcipi_numacc_x_potbig(p)) then
       rcipi_efoldmax = huge(1._kp)
       return
    endif
    
    xend = rcipi_x_endinf(p,alpha,beta)

    rcipi_efoldmax = rcipi_efold_primitive(xinimax,p,alpha,beta) &
         - rcipi_efold_primitive(xend,p,alpha,beta)

  end function rcipi_efoldmax
  
end module rcipisr
