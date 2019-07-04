!slow-roll function for string axion II inflation at increasing field
!values x > xVmax and x < xVmin, when V(xVmin)<= 0, such that inflation
!gracefully ends. The model exists only if both xVmax and xVmin exists
!and if V(xVmin)<= 0.
!
!!V(phi) = M^4 [1 - cos(x) + alpha x sin(x) + (1/2) alpha beta x^2]
!
!x=phi/mu
!
!
!
module saiii2sr
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use saiiicommon, only : saiii_norm_potential, saiii_norm_deriv_potential
  use saiiicommon, only : saiii_norm_deriv_second_potential
  use saiiicommon, only : saiii_epsilon_one, saiii_epsilon_two, saiii_epsilon_three
  use saiiicommon, only : saiii_x_potzero, saiii_x_potmax, saiii_efold_primitive
  use saiiicommon, only : saiii_x_epsoneunity, find_saiii_x_trajectory
  use saiiicommon, only : saiii_alpha_potneg, saiii_check_params, beta3, beta2
  
  implicit none

  
  private
  public saiii2_norm_potential, saiii2_norm_deriv_potential, saiii2_norm_deriv_second_potential
  public saiii2_epsilon_one, saiii2_epsilon_two, saiii2_epsilon_three
  public saiii2_x_endinf, saiii2_x_trajectory, saiii2_efold_primitive
  public saiii2_check_params
  public saiii2_numacc_xinimin, saiii2_numacc_efoldmax, saiii2_numacc_mumin

!numerical accuracy limitation
  real(kp), parameter :: epsnumacc = 10._kp*tolkp
  
  
contains


  function saiii2_check_params(alpha,beta,mu)
    implicit none
    logical :: saiii2_check_params
    real(kp), intent(in) :: alpha,beta,mu
    logical :: sane

    real(kp) :: alphabound

    sane = saiii_check_params(alpha,beta,mu) &
         .and. (beta.lt.beta2) &
         .and. (beta.gt.beta3)

    if (.not.sane) then

       saiii2_check_params = .false.

    else 

       alphabound = saiii_alpha_potneg(beta)

       saiii2_check_params = (abs(alpha).gt.abs(alphabound))

    endif               

  end function saiii2_check_params


  

  function saiii2_norm_potential(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii2_norm_potential
    real(kp), intent(in) :: x,alpha,beta,mu

    saiii2_norm_potential = saiii_norm_potential(x,alpha,beta,mu)
    
  end function saiii2_norm_potential


  
!derivative with respect to x (not phi!)  
  function saiii2_norm_deriv_potential(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii2_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta,mu


    saiii2_norm_deriv_potential =  saiii_norm_deriv_potential(x,alpha,beta,mu)

  end function saiii2_norm_deriv_potential

  

!second derivative with respect to x (not phi!)    
  function saiii2_norm_deriv_second_potential(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii2_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta,mu
    
    saiii2_norm_deriv_second_potential = saiii_norm_deriv_second_potential(x,alpha,beta,mu)

  end function saiii2_norm_deriv_second_potential



  
  function saiii2_epsilon_one(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii2_epsilon_one
    real(kp), intent(in) :: x,alpha,beta,mu
    
    saiii2_epsilon_one = saiii_epsilon_one(x,alpha,beta,mu)
    
  end function saiii2_epsilon_one
 
  

  
  function saiii2_epsilon_two(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii2_epsilon_two
    real(kp), intent(in) :: x,alpha,beta,mu
    
    saiii2_epsilon_two = saiii_epsilon_two(x,alpha,beta,mu)
    
  end function saiii2_epsilon_two



  
  function saiii2_epsilon_three(x,alpha,beta,mu)
    implicit none
    real(kp) :: saiii2_epsilon_three
    real(kp), intent(in) :: x,alpha,beta,mu
    
    saiii2_epsilon_three = saiii_epsilon_three(x,alpha,beta,mu)

  end function saiii2_epsilon_three



  
  function saiii2_x_endinf(alpha,beta,mu)
    implicit none
    real(kp) :: saiii2_x_endinf
    real(kp), intent(in) :: alpha,beta,mu

    real(kp), dimension(2) :: xepsone
    
    xepsone = saiii_x_epsoneunity(alpha,beta,mu)
    
    saiii2_x_endinf = xepsone(2)
    
  end function saiii2_x_endinf



!the closest possible to the top of the potential
  function saiii2_numacc_xinimin(alpha,beta,mu)
    implicit none
    real(kp) :: saiii2_numacc_xinimin
    real(kp), intent(in) :: alpha,beta,mu

    real(kp), parameter :: dx = 0.001_kp
    real(kp) :: dlnVodx
    
    real(kp), save :: xVmax = huge(1._kp)
    real(kp), save :: stoalpha = huge(1._kp)
    real(kp), save :: stobeta = huge(1._kp)
!$omp threadprivate(xVmax,stoalpha,stobeta)

    if (.not.saiii2_check_params(alpha,beta,mu)) then
       stop 'saiii2_numacc_xinimin: inflation at x > xVmax never ends!'
    endif
    
    if ((alpha.ne.stoalpha).or.(beta.ne.stobeta)) then
       xVmax = saiii_x_potmax(alpha,beta,mu)
       stoalpha = alpha
       stobeta = beta
    endif

    if (saiii2_norm_potential(xVmax+dx,alpha,beta,mu).lt.0._kp) then
       stop 'saiii2_numacc_xinimin: dx too large!'
    endif
    
    dlnVodx = abs(saiii2_norm_deriv_potential(xVmax+dx,alpha,beta,mu) &
         /saiii2_norm_potential(xVmax,alpha,beta,mu)/dx)

    saiii2_numacc_xinimin = xVmax + epsnumacc/dlnVodx
    
  end function saiii2_numacc_xinimin
    

!maximal number of efolds computable at current numerical accuracy
  function saiii2_numacc_efoldmax(alpha,beta,mu)
    implicit none
    real(kp) :: saiii2_numacc_efoldmax
    real(kp), intent(in) :: alpha,beta,mu

    real(kp) :: xend,xinimin
    
    if (.not.saiii2_check_params(alpha,beta,mu)) then
       stop 'saiii2_numacc_efoldmax: inflation at x > xVmax never ends!'
    endif
    
    xend = saiii2_x_endinf(alpha,beta,mu)
    xinimin = saiii2_numacc_xinimin(alpha,beta,mu)
    
    saiii2_numacc_efoldmax = -saiii2_efold_primitive(xend,alpha,beta,mu) &
         + saiii2_efold_primitive(xinimin,alpha,beta,mu)

    
  end function saiii2_numacc_efoldmax
 


!given alpha and beta what is the minimal value of mu to get efold inflation
!above numerical accuracy limit
  function saiii2_numacc_mumin(efold,alpha,beta)
    implicit none
    real(kp) :: saiii2_numacc_mumin
    real(kp), intent(in) :: efold,alpha,beta

    real(kp), parameter :: tolFind=tolkp
    type(transfert) :: saiii2Data

    real(kp), parameter :: mubig = 10000._kp
    real(kp), parameter :: musmall = 0.0001_kp
    
    real(kp) :: mini,maxi

    mini = musmall
    maxi = mubig
    
    saiii2Data%real1 = alpha
    saiii2Data%real2 = beta
    saiii2Data%real3 = efold

    saiii2_numacc_mumin = zbrent(find_saiii2_numacc_mumin,mini,maxi,tolFind,saiii2Data)

  end function saiii2_numacc_mumin



  function find_saiii2_numacc_mumin(mu,saiii2Data)
    implicit none
    real(kp) :: find_saiii2_numacc_mumin
    real(kp), intent(in) :: mu
    type(transfert), optional, intent(inout) :: saiii2Data

    real(kp) :: alpha, beta, Nwanted

    alpha = saiii2Data%real1
    beta = saiii2Data%real2
    Nwanted = saiii2Data%real3

    find_saiii2_numacc_mumin = saiii2_numacc_efoldmax(alpha,beta,mu) - Nwanted
    
  end function find_saiii2_numacc_mumin

  
  
  function saiii2_efold_primitive(x,alpha,beta,mu)
    implicit none
    real(kp), intent(in) :: x, alpha,beta,mu
    real(kp) :: saiii2_efold_primitive

    if (x.lt.saiii2_numacc_xinimin(alpha,beta,mu)) then
       write(*,*)'saiii2_efold_primitive: x-xVmax too small!'
       write(*,*)'x= alpha= mu= ',x,alpha,beta,mu
       stop
    endif
    
    saiii2_efold_primitive = saiii_efold_primitive(x,alpha,beta,mu)

  end function saiii2_efold_primitive


  
!  returns x at bfold=-efolds before the end of inflation
  function saiii2_x_trajectory(bfold,xend,alpha,beta,mu)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,beta,mu
    real(kp) :: saiii2_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: saiiiData

    real(kp) :: xinimin
    
    xinimin = saiii2_numacc_xinimin(alpha,beta,mu)

    saiiiData%real1 = alpha
    saiiiData%real2 = beta
    saiiiData%real3 = mu
    saiiiData%real4 = -bfold + saiii2_efold_primitive(xend,alpha,beta,mu)

    saiii2_x_trajectory = zbrent(find_saiii_x_trajectory,xinimin,xend,tolFind,saiiiData)

  end function saiii2_x_trajectory

  
end module saiii2sr
