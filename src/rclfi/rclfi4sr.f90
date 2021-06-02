!Radiatively corrected large field inflation when the potential is a
!growing monotonous function of the field everywhere. The potential
!has no extremum.
!
!
!V(phi) = M^4 [x^p + alpha x^4 ln(x) ]
!
!x = phi/mu
!
!with:
!  p>4 and -p(p-4)/4exp(2-p/4) < alpha < 0
!OR
!  p<4 and 0 < alpha < -p(p-4)/4exp(2-p/4)
!

module rclfi4sr
  use infprec, only : kp,tolkp,toldp,transfert
  use inftools, only : zbrent, easydverk
  
  use rclficommon, only : rclfi_check_params
  use rclficommon, only : rclfi_norm_potential, rclfi_norm_deriv_potential
  use rclficommon, only : rclfi_norm_deriv_second_potential, rclfi_epsilon_one
  use rclficommon, only : rclfi_epsilon_two, rclfi_epsilon_three
  use rclficommon, only : find_rclfi_x_endinf, rclfi_alpha_one
  use rclficommon, only : find_rclfi_efold_primitive, rclfi_x_derivpotzero

  implicit none

  private

  

!essentially LFI at large field value with p>=4, above this efold_primitive
!becomes smaller than numerical accuracy
  real(kp), parameter :: rclfixbig = 10000._kp
  
  public rclfixBig
  public rclfi4_check_params, rclfi4_norm_deriv_second_potential
  public rclfi4_norm_potential, rclfi4_norm_deriv_potential
  public rclfi4_epsilon_one, rclfi4_epsilon_two, rclfi4_epsilon_three
  public rclfi4_x_endinf, rclfi4_efold_primitive, rclfi4_x_trajectory

  
contains


  
  function rclfi4_check_params(p,alpha,mu)
    logical :: rclfi4_check_params
    real(kp), intent(in) :: p,alpha,mu

    real(kp) :: alphadVo

    alphadVo = rclfi_alpha_one(p)
    
    rclfi4_check_params = ((p.gt.4._kp).and.(alpha.le.0._kp).and.(alpha.ge.alphadVo)) &
         .or.((p.lt.4._kp).and.(alpha.ge.0._kp).and.(alpha.le.alphadVo))
    
    rclfi4_check_params = rclfi4_check_params.and.rclfi_check_params(p,alpha,mu)
    
  end function rclfi4_check_params
  

  
  function rclfi4_norm_potential(x,p,alpha,mu)
    implicit none
    real(kp) :: rclfi4_norm_potential
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi4_norm_potential = rclfi_norm_potential(x,p,alpha,mu)

  end function rclfi4_norm_potential




  function rclfi4_norm_deriv_potential(x,p,alpha,mu)
    implicit none
    real(kp) :: rclfi4_norm_deriv_potential
    real(kp), intent(in) :: x,p,alpha,mu

   rclfi4_norm_deriv_potential = rclfi_norm_deriv_potential(x,p,alpha,mu)

  end function rclfi4_norm_deriv_potential




  function rclfi4_norm_deriv_second_potential(x,p,alpha,mu)
    implicit none
    real(kp) :: rclfi4_norm_deriv_second_potential
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi4_norm_deriv_second_potential = rclfi_norm_deriv_second_potential(x,p,alpha,mu)

  end function rclfi4_norm_deriv_second_potential




  function rclfi4_epsilon_one(x,p,alpha,mu)    
    implicit none
    real(kp) :: rclfi4_epsilon_one
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi4_epsilon_one = rclfi_epsilon_one(x,p,alpha,mu)
    
  end function rclfi4_epsilon_one



  function rclfi4_epsilon_two(x,p,alpha,mu)    
    implicit none
    real(kp) :: rclfi4_epsilon_two
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi4_epsilon_two = rclfi_epsilon_two(x,p,alpha,mu)
    
  end function rclfi4_epsilon_two



  function rclfi4_epsilon_three(x,p,alpha,mu)    
    implicit none
    real(kp) :: rclfi4_epsilon_three
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi4_epsilon_three = rclfi_epsilon_three(x,p,alpha,mu)
    
  end function rclfi4_epsilon_three


 
  

!returns x at the end of inflation defined as epsilon1=1
  function rclfi4_x_endinf(p,alpha,mu)
    implicit none
    real(kp), intent(in) :: p,alpha,mu
    real(kp) :: rclfi4_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rclfiData

    mini = epsilon(1._kp)
    maxi = rclfixBig

    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    
    rclfi4_x_endinf = zbrent(find_rclfi_x_endinf,mini,maxi,tolFind,rclfiData)
    
  end function rclfi4_x_endinf



!this is integral[V(phi)/V'(phi) dphi]
  function rclfi4_efold_primitive(x,p,alpha,mu)
    implicit none
    real(kp), intent(in) :: x,p,alpha,mu
    real(kp) :: rclfi4_efold_primitive

    type(transfert) :: rclfiData
!too long to integrate in QUADPREC
    real(kp), parameter :: tolInt = 10._kp*max(tolkp,toldp)
    integer, parameter :: neq = 1

    real(kp) :: xvar
    real(kp), dimension(neq) :: yvar

  
    xvar = 0._kp
    yvar(1) = 0._kp
    
    rclfiData%real1 = alpha
    rclfiData%real2 = p
    
    call easydverk(neq,find_rclfi_efold_primitive,xvar,yvar,x,tolInt,rclfiData)

    rclfi4_efold_primitive = yvar(1) * mu * mu
    
  end function rclfi4_efold_primitive


  
  
!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function rclfi4_x_trajectory(bfold,xend,p,alpha,mu)
    implicit none
    real(kp), intent(in) :: bfold, alpha, p, mu, xend
    real(kp) :: rclfi4_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rclfiData


    mini = xend
    maxi = rclfixBig

  
    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    rclfiData%real4 = -bfold + rclfi4_efold_primitive(xend,p,alpha,mu)

    rclfi4_x_trajectory = zbrent(find_rclfi4_x_trajectory,mini,maxi,tolFind,rclfiData)
    
  end function rclfi4_x_trajectory


  
  function find_rclfi4_x_trajectory(x,rclfiData)
    implicit none
    real(kp) :: find_rclfi4_x_trajectory
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: alpha, p, mu, NplusNuend

    alpha = rclfiData%real1
    p = rclfiData%real2
    mu = rclfiData%real3    
    NplusNuend = rclfiData%real4

    find_rclfi4_x_trajectory = rclfi4_efold_primitive(x,p,alpha,mu) - NplusNuend

  end function find_rclfi4_x_trajectory

  
end module rclfi4sr
