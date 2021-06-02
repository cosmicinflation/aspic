!Radiatively corrected large field inflation in the large field
!region; at field value larger than the minimum of the potential and
!where it is positive. Parameters are such that the potential has a
!minimum, if not, see rclfi4.
!
!V(phi) = M^4 [x^p + alpha x^4 ln(x) ]
!
!x = phi/mu
!
!with:
!  p>4 and alpha>0
!OR
!  p<4 and alpha<-e(p-4) [<0]
!OR
!  p<4 and alpha> -e(p-4) [>-p(p-4)/4 exp(2-p/4)>0]
!

module rclfi3sr
  use infprec, only : kp,tolkp,toldp,transfert
  use inftools, only : zbrent, easydverk
  
  use rclficommon, only : rclfi_check_params
  use rclficommon, only : rclfi_norm_potential, rclfi_norm_deriv_potential
  use rclficommon, only : rclfi_norm_deriv_second_potential, rclfi_epsilon_one
  use rclficommon, only : rclfi_epsilon_two, rclfi_epsilon_three, rclfi_x_potzero
  use rclficommon, only : find_rclfi_x_endinf, rclfi_alpha_zero, rclfi_alpha_one
  use rclficommon, only : find_rclfi_efold_primitive, rclfi_x_derivpotzero

  implicit none

  private

  
!numerical accuracy limitation
  real(kp), parameter :: epsnumacc = 100._kp*tolkp

  !essentially LFI at large field value with p>=4, above this eps1
!becomes smaller than machine accuracy
  real(kp), parameter :: xpotzeroMax = 500._kp
  real(kp), parameter :: rclfixbig = 100._kp*xpotzeromax
  
  public rclfixBig
  public rclfi3_check_params, rclfi3_norm_deriv_second_potential
  public rclfi3_norm_potential, rclfi3_norm_deriv_potential
  public rclfi3_epsilon_one, rclfi3_epsilon_two, rclfi3_epsilon_three
  public rclfi3_x_endinf, rclfi3_efold_primitive, rclfi3_x_trajectory
  public rclfi3_x_potzero, rclfi3_numacc_alphamin, rclfi3_numacc_pmin

  
contains


  
  function rclfi3_check_params(p,alpha,mu)
    logical :: rclfi3_check_params
    real(kp), intent(in) :: p,alpha,mu

    real(kp) :: alphaVo

    alphaVo = rclfi_alpha_zero(p)
    
    rclfi3_check_params = ((p.gt.4._kp).and.(alpha.ge.0._kp)) &
         .or. ((p.gt.4._kp).and.(alpha.le.alphaVo)) &
         .or. ((p.lt.4._kp).and.(alpha.ge.alphaVo))
    
    rclfi3_check_params = rclfi3_check_params.and.rclfi_check_params(p,alpha,mu)
    
  end function rclfi3_check_params
  

  
  function rclfi3_norm_potential(x,p,alpha,mu)
    implicit none
    real(kp) :: rclfi3_norm_potential
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi3_norm_potential = rclfi_norm_potential(x,p,alpha,mu)

  end function rclfi3_norm_potential




  function rclfi3_norm_deriv_potential(x,p,alpha,mu)
    implicit none
    real(kp) :: rclfi3_norm_deriv_potential
    real(kp), intent(in) :: x,p,alpha,mu

   rclfi3_norm_deriv_potential = rclfi_norm_deriv_potential(x,p,alpha,mu)

  end function rclfi3_norm_deriv_potential




  function rclfi3_norm_deriv_second_potential(x,p,alpha,mu)
    implicit none
    real(kp) :: rclfi3_norm_deriv_second_potential
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi3_norm_deriv_second_potential = rclfi_norm_deriv_second_potential(x,p,alpha,mu)

  end function rclfi3_norm_deriv_second_potential




  function rclfi3_epsilon_one(x,p,alpha,mu)    
    implicit none
    real(kp) :: rclfi3_epsilon_one
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi3_epsilon_one = rclfi_epsilon_one(x,p,alpha,mu)
    
  end function rclfi3_epsilon_one



  function rclfi3_epsilon_two(x,p,alpha,mu)    
    implicit none
    real(kp) :: rclfi3_epsilon_two
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi3_epsilon_two = rclfi_epsilon_two(x,p,alpha,mu)
    
  end function rclfi3_epsilon_two



  function rclfi3_epsilon_three(x,p,alpha,mu)    
    implicit none
    real(kp) :: rclfi3_epsilon_three
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi3_epsilon_three = rclfi_epsilon_three(x,p,alpha,mu)
    
  end function rclfi3_epsilon_three


 
!non vanishing x at which the potential becomes negative or null
  function rclfi3_x_potzero(p,alpha,mu)
    implicit none
    real(kp) :: rclfi3_x_potzero
    real(kp), intent(in) :: p,alpha,mu

    real(kp), dimension(2) :: xVzero

    if (.not.rclfi3_check_params(p,alpha,mu)) then
       stop 'rclfi3_x_potzero: no zeros!'
    endif
    
    xVzero = rclfi_x_potzero(p,alpha,mu)

    rclfi3_x_potzero = xVzero(2)
    
  end function rclfi3_x_potzero

 
  

!returns x at the end of inflation defined as epsilon1=1
  function rclfi3_x_endinf(p,alpha,mu)
    implicit none
    real(kp), intent(in) :: p,alpha,mu
    real(kp) :: rclfi3_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rclfiData

    mini = rclfi3_x_potzero(p,alpha,mu)+ epsilon(1._kp)
    maxi = rclfixBig

    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    
    rclfi3_x_endinf = zbrent(find_rclfi_x_endinf,mini,maxi,tolFind,rclfiData)
    
  end function rclfi3_x_endinf



!returns alphamin<0 such that for alpha > alphamin, we have xpotzero < xbig
  function rclfi3_numacc_alphamin_x(xbig,p)
    implicit none
    real(kp) :: rclfi3_numacc_alphamin_x
    real(kp), intent(in) :: xbig,p

    if (p.lt.4._kp) then
       stop 'rclfi3_numacc_alphamin_x: p<4!'
    endif
    
    rclfi3_numacc_alphamin_x = -(p-4._kp)*xbig**(p-4._kp)/log(xbig**(p-4._kp))

    if (rclfi3_numacc_alphamin_x.gt.rclfi_alpha_zero(p)) then
       write(*,*)'alphamin= alphaVo= ',rclfi3_numacc_alphamin_x,rclfi_alpha_zero(p)
       stop 'p too close to 4'
    endif
    
  end function rclfi3_numacc_alphamin_x



  function rclfi3_numacc_alphamin(p)
    implicit none
    real(kp) :: rclfi3_numacc_alphamin
    real(kp), intent(in) :: p

    rclfi3_numacc_alphamin = rclfi3_numacc_alphamin_x(xpotzeromax,p)

  end function rclfi3_numacc_alphamin


  
!returns pmin>4 such that for all values of alpha>alphaVo we have
!xpotzero < xbig
  function rclfi3_numacc_pmin_x(xbig)
    implicit none
    real(kp) :: rclfi3_numacc_pmin_x
    real(kp), intent(in) :: xbig

    rclfi3_numacc_pmin_x = 4._kp + 1._kp/log(xbig)

  end function rclfi3_numacc_pmin_x


  function rclfi3_numacc_pmin()
    implicit none
    real(kp) :: rclfi3_numacc_pmin

    rclfi3_numacc_pmin = rclfi3_numacc_pmin_x(xpotzeromax)
    
  end function rclfi3_numacc_pmin
  

!this is integral[V(phi)/V'(phi) dphi]
  function rclfi3_efold_primitive(x,p,alpha,mu)
    implicit none
    real(kp), intent(in) :: x,p,alpha,mu
    real(kp) :: rclfi3_efold_primitive

    type(transfert) :: rclfiData
!too long to integrate in QUADPREC
    real(kp), parameter :: tolInt = max(tolkp,toldp)
    integer, parameter :: neq = 1

    real(kp) :: xpotzero
    real(kp) :: xvar
    real(kp), dimension(neq) :: yvar

    
    xpotzero = rclfi3_x_potzero(p,alpha,mu)
    
    if (x.le.xpotzero) then
       stop 'rclfi3_efold_primitive: x not in ]xpotzero,+oo[!'
    endif

    xvar = xpotzero
    yvar(1) = 0._kp
    
    rclfiData%real1 = alpha
    rclfiData%real2 = p
    
    call easydverk(neq,find_rclfi_efold_primitive,xvar,yvar,x,tolInt,rclfiData)

    rclfi3_efold_primitive = yvar(1) * mu * mu
    
  end function rclfi3_efold_primitive


  
  
!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function rclfi3_x_trajectory(bfold,xend,p,alpha,mu)
    implicit none
    real(kp), intent(in) :: bfold, alpha, p, mu, xend
    real(kp) :: rclfi3_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rclfiData


    mini = xend
    maxi = rclfixBig

  
    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    rclfiData%real4 = -bfold + rclfi3_efold_primitive(xend,p,alpha,mu)

    rclfi3_x_trajectory = zbrent(find_rclfi3_x_trajectory,mini,maxi,tolFind,rclfiData)
    
  end function rclfi3_x_trajectory


  
  function find_rclfi3_x_trajectory(x,rclfiData)
    implicit none
    real(kp) :: find_rclfi3_x_trajectory
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: alpha, p, mu, NplusNuend

    alpha = rclfiData%real1
    p = rclfiData%real2
    mu = rclfiData%real3    
    NplusNuend = rclfiData%real4

    find_rclfi3_x_trajectory = rclfi3_efold_primitive(x,p,alpha,mu) - NplusNuend

  end function find_rclfi3_x_trajectory

  
end module rclfi3sr
