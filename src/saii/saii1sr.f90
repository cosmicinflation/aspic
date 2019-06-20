!slow-roll function for string axion I inflation at decreasing field
!values x < xVmax
!
!V(phi) = M^4 [1 - cos(x) + alpha x sin(x)]
!
!x=phi/mu
!
!
!
module saii1sr
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use saiicommon, only : saii_norm_potential, saii_norm_deriv_potential
  use saiicommon, only : saii_norm_deriv_second_potential
  use saiicommon, only : saii_epsilon_one, saii_epsilon_two, saii_epsilon_three
  use saiicommon, only : saii_x_potzero, saii_x_potmax, saii_efold_primitive
  use saiicommon, only : saii_x_epsoneunity, find_saii_x_trajectory

  
  implicit none

  
  private
  public saii1_norm_potential, saii1_norm_deriv_potential, saii1_norm_deriv_second_potential
  public saii1_epsilon_one, saii1_epsilon_two, saii1_epsilon_three
  public saii1_x_endinf, saii1_x_trajectory, saii1_efold_primitive
  public saii1_numacc_xinimax, saii1_numacc_efoldmax, saii1_numacc_mumin

!numerical accuracy limitation
  real(kp), parameter :: epsnumacc = 10._kp*tolkp
  
  
contains


  function saii1_norm_potential(x,alpha,mu)
    implicit none
    real(kp) :: saii1_norm_potential
    real(kp), intent(in) :: x,alpha,mu

    saii1_norm_potential = saii_norm_potential(x,alpha,mu)
    
  end function saii1_norm_potential


  
!derivative with respect to x (not phi!)  
  function saii1_norm_deriv_potential(x,alpha,mu)
    implicit none
    real(kp) :: saii1_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,mu


    saii1_norm_deriv_potential =  saii_norm_deriv_potential(x,alpha,mu)

  end function saii1_norm_deriv_potential

  

!second derivative with respect to x (not phi!)    
  function saii1_norm_deriv_second_potential(x,alpha,mu)
    implicit none
    real(kp) :: saii1_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,mu
    
    saii1_norm_deriv_second_potential = saii_norm_deriv_second_potential(x,alpha,mu)

  end function saii1_norm_deriv_second_potential



  
  function saii1_epsilon_one(x,alpha,mu)
    implicit none
    real(kp) :: saii1_epsilon_one
    real(kp), intent(in) :: x,alpha,mu
    
    saii1_epsilon_one = saii_epsilon_one(x,alpha,mu)
    
  end function saii1_epsilon_one
 
  

  
  function saii1_epsilon_two(x,alpha,mu)
    implicit none
    real(kp) :: saii1_epsilon_two
    real(kp), intent(in) :: x,alpha,mu
    
    saii1_epsilon_two = saii_epsilon_two(x,alpha,mu)
    
  end function saii1_epsilon_two



  
  function saii1_epsilon_three(x,alpha,mu)
    implicit none
    real(kp) :: saii1_epsilon_three
    real(kp), intent(in) :: x,alpha,mu
    
    saii1_epsilon_three = saii_epsilon_three(x,alpha,mu)

  end function saii1_epsilon_three



  
  function saii1_x_endinf(alpha,mu)
    implicit none
    real(kp) :: saii1_x_endinf
    real(kp), intent(in) :: alpha, mu

    real(kp), dimension(2) :: xepsone
    
    xepsone = saii_x_epsoneunity(alpha,mu)
    
    saii1_x_endinf = xepsone(1)
    
  end function saii1_x_endinf



!the closest possible to the top of the potential
  function saii1_numacc_xinimax(alpha,mu)
    implicit none
    real(kp) :: saii1_numacc_xinimax
    real(kp), intent(in) :: alpha,mu

    real(kp), parameter :: dx = 0.001_kp
    real(kp) :: dlnVodx
    
    real(kp), save :: xVmax = huge(1._kp)
    real(kp), save :: stoalpha = huge(1._kp)
!$omp threadprivate(xVmax,stoalpha)

    if (alpha.ne.stoalpha) then
       xVmax = saii_x_potmax(alpha,mu)
       stoalpha = alpha
    endif

    if (saii1_norm_potential(xVmax-dx,alpha,mu).lt.0._kp) then
       stop 'saii1_numacc_xinimax: dx too large!'
    endif
    
    dlnVodx = abs(saii1_norm_deriv_potential(xVmax-dx,alpha,mu) &
         /saii1_norm_potential(xVmax,alpha,mu)/dx)

    saii1_numacc_xinimax = xVmax - epsnumacc/dlnVodx
    
  end function saii1_numacc_xinimax
    

!maximal number of efolds computable at current numerical accuracy
  function saii1_numacc_efoldmax(alpha,mu)
    implicit none
    real(kp) :: saii1_numacc_efoldmax
    real(kp), intent(in) :: alpha,mu

    real(kp) :: xend,xinimax
    
    
    xend = saii1_x_endinf(alpha,mu)
    xinimax = saii1_numacc_xinimax(alpha,mu)
    
    saii1_numacc_efoldmax = -saii1_efold_primitive(xend,alpha,mu) &
         + saii1_efold_primitive(xinimax,alpha,mu)

    
  end function saii1_numacc_efoldmax
 


!given alpha, what is the minimal value of mu to get efold inflation
!above numerical accuracy limit
  function saii1_numacc_mumin(efold,alpha)
    implicit none
    real(kp) :: saii1_numacc_mumin
    real(kp), intent(in) :: efold,alpha

    real(kp), parameter :: tolFind=tolkp
    type(transfert) :: saii1Data

    real(kp), parameter :: mubig = 10000._kp
    real(kp), parameter :: musmall = 0.0001_kp
    
    real(kp) :: mini,maxi

    mini = musmall
    maxi = mubig
    
    saii1Data%real1 = alpha
    saii1Data%real2 = efold

    saii1_numacc_mumin = zbrent(find_saii1_numacc_mumin,mini,maxi,tolFind,saii1Data)

  end function saii1_numacc_mumin



  function find_saii1_numacc_mumin(mu,saii1Data)
    implicit none
    real(kp) :: find_saii1_numacc_mumin
    real(kp), intent(in) :: mu
    type(transfert), optional, intent(inout) :: saii1Data

    real(kp) :: alpha, Nwanted

    alpha = saii1Data%real1
    Nwanted = saii1Data%real2

    find_saii1_numacc_mumin = saii1_numacc_efoldmax(alpha,mu) - Nwanted
    
  end function find_saii1_numacc_mumin

  
  
  function saii1_efold_primitive(x,alpha,mu)
    implicit none
    real(kp), intent(in) :: x, alpha,mu
    real(kp) :: saii1_efold_primitive

    if (x.gt.saii1_numacc_xinimax(alpha,mu)) then
       write(*,*)'saii1_efold_primitive: xVmax-x too small!'
       write(*,*)'x= alpha= mu= ',x,alpha,mu
       stop
    endif
    
    saii1_efold_primitive = saii_efold_primitive(x,alpha,mu)

  end function saii1_efold_primitive


  
!  returns x at bfold=-efolds before the end of inflation
  function saii1_x_trajectory(bfold,xend,alpha,mu)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,mu
    real(kp) :: saii1_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: saiiData

    real(kp) :: xinimax
    
    xinimax = saii1_numacc_xinimax(alpha,mu)

    saiiData%real1 = alpha
    saiiData%real2 = mu
    saiiData%real3 = -bfold + saii1_efold_primitive(xend,alpha,mu)

    saii1_x_trajectory = zbrent(find_saii_x_trajectory,xend,xinimax,tolFind,saiiData)

  end function saii1_x_trajectory

  
end module saii1sr
