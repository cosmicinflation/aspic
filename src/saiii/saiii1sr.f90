!slow-roll function for string axion II inflation at decreasing field
!values x < xVmax, when xVmax exists
!
!V(phi) = M^4 [1 - cos(x) + alpha x sin(x) + (1/2) alpha beta x^2]
!
!x=phi/mu
!
!
!
module saiii1sr
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use saiiicommon, only : saiii_norm_potential, saiii_norm_deriv_potential
  use saiiicommon, only : saiii_norm_deriv_second_potential
  use saiiicommon, only : saiii_epsilon_one, saiii_epsilon_two, saiii_epsilon_three
  use saiiicommon, only : saiii_x_potzero, saiii_x_potmax, saiii_efold_primitive
  use saiiicommon, only : saiii_x_epsoneunity, find_saiii_x_trajectory
  use saiiicommon, only : saiii_check_minima

  
  implicit none

  
  private
  public saiii1_norm_potential, saiii1_norm_deriv_potential, saiii1_norm_deriv_second_potential
  public saiii1_epsilon_one, saiii1_epsilon_two, saiii1_epsilon_three
  public saiii1_x_endinf, saiii1_x_trajectory, saiii1_efold_primitive
  public saiii1_numacc_xinimax, saiii1_numacc_efoldmax, saiii1_numacc_mumin

!numerical accuracy limitation
  real(kp), parameter :: epsnumacc = 10._kp*tolkp
  
  
contains


  function saiii1_check_params(alpha,beta,mu)
    implicit none
    logical :: saiii1_check_params
    real(kp), intent(in) :: alpha,beta,mu

    saiii1_check_params = saiii_check_minima(alpha,beta,mu)
    
  end function saiii1_check_params

  
  function saiii1_norm_potential(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii1_norm_potential
    real(kp), intent(in) :: x,alpha,beta,mu

    saiii1_norm_potential = saiii_norm_potential(x,alpha,beta,mu)
    
  end function saiii1_norm_potential


  
!derivative with respect to x (not phi!)  
  function saiii1_norm_deriv_potential(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii1_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta,mu

    saiii1_norm_deriv_potential =  saiii_norm_deriv_potential(x,alpha,beta,mu)

  end function saiii1_norm_deriv_potential

  

!second derivative with respect to x (not phi!)    
  function saiii1_norm_deriv_second_potential(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii1_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta,mu
    
    saiii1_norm_deriv_second_potential = saiii_norm_deriv_second_potential(x,alpha,beta,mu)

  end function saiii1_norm_deriv_second_potential



  
  function saiii1_epsilon_one(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii1_epsilon_one
    real(kp), intent(in) :: x,alpha,beta,mu
    
    saiii1_epsilon_one = saiii_epsilon_one(x,alpha,beta,mu)
    
  end function saiii1_epsilon_one
 
  

  
  function saiii1_epsilon_two(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii1_epsilon_two
    real(kp), intent(in) :: x,alpha,beta,mu
    
    saiii1_epsilon_two = saiii_epsilon_two(x,alpha,beta,mu)
    
  end function saiii1_epsilon_two



  
  function saiii1_epsilon_three(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii1_epsilon_three
    real(kp), intent(in) :: x,alpha,beta,mu
    
    saiii1_epsilon_three = saiii_epsilon_three(x,alpha,beta,mu)

  end function saiii1_epsilon_three



  
  function saiii1_x_endinf(alpha,beta,mu)
    implicit none
    real(kp) :: saiii1_x_endinf
    real(kp), intent(in) :: alpha, beta, mu

    real(kp), dimension(2) :: xepsone
    
    xepsone = saiii_x_epsoneunity(alpha,beta,mu)
    
    saiii1_x_endinf = xepsone(1)
    
  end function saiii1_x_endinf



!the closest possible to the top of the potential
  function saiii1_numacc_xinimax(alpha,beta,mu)
    implicit none
    real(kp) :: saiii1_numacc_xinimax
    real(kp), intent(in) :: alpha,beta,mu

    real(kp), parameter :: dx = 0.001_kp
    real(kp) :: dlnVodx
    
    real(kp), save :: xVmax = huge(1._kp)
    real(kp), save :: stoalpha = huge(1._kp)
    real(kp), save :: stobeta = huge(1._kp)
!$omp threadprivate(xVmax,stoalpha,stobeta)


    if (.not.saiii1_check_params(alpha,beta,mu)) then
       stop 'saiii1_numacc_xinimax: potential has no maximum!'
    endif
    
    
    if ((alpha.ne.stoalpha).or.(beta.ne.stobeta)) then
       xVmax = saiii_x_potmax(alpha,beta,mu)
       stoalpha = alpha
       stobeta = beta
    endif

    if (saiii1_norm_potential(xVmax-dx,alpha,beta,mu).lt.0._kp) then
       write(*,*)'alpha= beta= ',alpha,beta
       stop 'saiii1_numacc_xinimax: dx too large!'
    endif
    
    dlnVodx = abs(saiii1_norm_deriv_potential(xVmax-dx,alpha,beta,mu) &
         /saiii1_norm_potential(xVmax,alpha,beta,mu)/dx)

    saiii1_numacc_xinimax = xVmax - epsnumacc/dlnVodx
    
  end function saiii1_numacc_xinimax
    

!maximal number of efolds computable at current numerical accuracy
  function saiii1_numacc_efoldmax(alpha,beta,mu)
    implicit none
    real(kp) :: saiii1_numacc_efoldmax
    real(kp), intent(in) :: alpha,beta,mu

    real(kp) :: xend,xinimax

    if (.not.saiii1_check_params(alpha,beta,mu)) then
       write(*,*)'alpha= beta= ',alpha,beta
       stop 'saiii1_numacc_efoldmax: potential has no maximum!'
    endif

    xend = saiii1_x_endinf(alpha,beta,mu)

    xinimax = saiii1_numacc_xinimax(alpha,beta,mu)
    
    saiii1_numacc_efoldmax = -saiii1_efold_primitive(xend,alpha,beta,mu) &
         + saiii1_efold_primitive(xinimax,alpha,beta,mu)

    
  end function saiii1_numacc_efoldmax
 


!given alpha, beta what is the minimal value of mu to get efold inflation
!above numerical accuracy limit
  function saiii1_numacc_mumin(efold,alpha,beta)
    implicit none
    real(kp) :: saiii1_numacc_mumin
    real(kp), intent(in) :: efold,alpha,beta

    real(kp), parameter :: tolFind=tolkp
    type(transfert) :: saiii1Data

    real(kp), parameter :: mubig = 10000._kp
    real(kp), parameter :: musmall = 0.0001_kp
    
    real(kp) :: mini,maxi

    mini = musmall
    maxi = mubig
    
    saiii1Data%real1 = alpha
    saiii1Data%real2 = beta
    saiii1Data%real3 = efold

    saiii1_numacc_mumin = zbrent(find_saiii1_numacc_mumin,mini,maxi,tolFind,saiii1Data)

  end function saiii1_numacc_mumin



  function find_saiii1_numacc_mumin(mu,saiii1Data)
    implicit none
    real(kp) :: find_saiii1_numacc_mumin
    real(kp), intent(in) :: mu
    type(transfert), optional, intent(inout) :: saiii1Data

    real(kp) :: alpha, beta, Nwanted

    alpha = saiii1Data%real1
    beta = saiii1Data%real2
    Nwanted = saiii1Data%real3

    find_saiii1_numacc_mumin = saiii1_numacc_efoldmax(alpha,beta,mu) - Nwanted
    
  end function find_saiii1_numacc_mumin

  
  
  function saiii1_efold_primitive(x,alpha,beta,mu)
    implicit none
    real(kp), intent(in) :: x, alpha,beta,mu
    real(kp) :: saiii1_efold_primitive

    if (x.gt.saiii1_numacc_xinimax(alpha,beta,mu)) then
       write(*,*)'saiii1_efold_primitive: xVmax-x too small!'
       write(*,*)'x= alpha= mu= ',x,alpha,beta,mu
       stop
    endif
    
    saiii1_efold_primitive = saiii_efold_primitive(x,alpha,beta,mu)

  end function saiii1_efold_primitive


  
!  returns x at bfold=-efolds before the end of inflation
  function saiii1_x_trajectory(bfold,xend,alpha,beta,mu)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,beta,mu
    real(kp) :: saiii1_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: saiiiData

    real(kp) :: xinimax
    
    xinimax = saiii1_numacc_xinimax(alpha,beta,mu)

    saiiiData%real1 = alpha
    saiiiData%real2 = beta
    saiiiData%real3 = mu
    saiiiData%real4 = -bfold + saiii1_efold_primitive(xend,alpha,beta,mu)

    saiii1_x_trajectory = zbrent(find_saiii_x_trajectory,xend,xinimax,tolFind,saiiiData)

  end function saiii1_x_trajectory

  
end module saiii1sr
