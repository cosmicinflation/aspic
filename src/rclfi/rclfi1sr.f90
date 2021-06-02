!Radiatively corrected large field inflation in the small
!field region (before the first maximum of the potential)
!
!V(phi) = M^4 [x^p + alpha x^4 ln(x) ]
!
!x = phi/mu
!
!with:
!  p>4 and alpha<-p(p-4)/4 exp(2-p/4)
!OR
!  p<4 and alpha<0
!OR
!  p<4 and alpha>-p(p-4)/4 exp(2-p/4)
!

module rclfi1sr
  use infprec, only : kp,tolkp,toldp,transfert
  use inftools, only : zbrent, easydverk
  
  use rclficommon, only : rclfi_check_params
  use rclficommon, only : rclfi_norm_potential, rclfi_norm_deriv_potential
  use rclficommon, only : rclfi_norm_deriv_second_potential, rclfi_epsilon_one
  use rclficommon, only : rclfi_epsilon_two, rclfi_epsilon_three
  use rclficommon, only : find_rclfi_x_endinf, rclfi_alpha_zero, rclfi_alpha_one
  use rclficommon, only : find_rclfi_efold_primitive, rclfi_x_derivpotzero

  implicit none

  private

  public rclfi1_check_params, rclfi1_norm_deriv_second_potential
  public rclfi1_norm_potential, rclfi1_norm_deriv_potential
  public rclfi1_epsilon_one, rclfi1_epsilon_two, rclfi1_epsilon_three
  public rclfi1_x_endinf, rclfi1_efold_primitive, rclfi1_x_trajectory
  public rclfi1_x_potmax, rclfi1_max_xpotmax, rclfi1_numacc_alphamax
  public rclfi1_numacc_xinimax, rclfi1_numacc_efoldmax, rclfi1_numacc_mumin
  public rclfi1_numacc_pmax

!numerical accuracy limitation
  real(kp), parameter :: epsnumacc = 10._kp*tolkp

!we do not consider inflation occuring at x< xpotmax < xpotmaxmin
  real(kp), parameter :: xpotmaxmin = 0.001_kp
  
contains


  
  function rclfi1_check_params(p,alpha,mu)
    logical :: rclfi1_check_params
    real(kp), intent(in) :: p,alpha,mu

    real(kp) :: alphadV

    alphadV = rclfi_alpha_one(p)
    
    rclfi1_check_params = ((p.gt.4._kp).and.(alpha.le.alphadV)) &
         .or. ((p.lt.4._kp).and.(alpha.lt.0._kp)) &
         .or. ((p.lt.4._kp).and.(alpha.ge.alphadV))
    
    rclfi1_check_params = rclfi1_check_params.and.rclfi_check_params(p,alpha,mu)
    
  end function rclfi1_check_params
  

  
  function rclfi1_norm_potential(x,p,alpha,mu)
    implicit none
    real(kp) :: rclfi1_norm_potential
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi1_norm_potential = rclfi_norm_potential(x,p,alpha,mu)

  end function rclfi1_norm_potential




  function rclfi1_norm_deriv_potential(x,p,alpha,mu)
    implicit none
    real(kp) :: rclfi1_norm_deriv_potential
    real(kp), intent(in) :: x,p,alpha,mu

   rclfi1_norm_deriv_potential = rclfi_norm_deriv_potential(x,p,alpha,mu)

  end function rclfi1_norm_deriv_potential




  function rclfi1_norm_deriv_second_potential(x,p,alpha,mu)
    implicit none
    real(kp) :: rclfi1_norm_deriv_second_potential
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi1_norm_deriv_second_potential = rclfi_norm_deriv_second_potential(x,p,alpha,mu)

  end function rclfi1_norm_deriv_second_potential




  function rclfi1_epsilon_one(x,p,alpha,mu)    
    implicit none
    real(kp) :: rclfi1_epsilon_one
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi1_epsilon_one = rclfi_epsilon_one(x,p,alpha,mu)
    
  end function rclfi1_epsilon_one



  function rclfi1_epsilon_two(x,p,alpha,mu)    
    implicit none
    real(kp) :: rclfi1_epsilon_two
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi1_epsilon_two = rclfi_epsilon_two(x,p,alpha,mu)
    
  end function rclfi1_epsilon_two



  function rclfi1_epsilon_three(x,p,alpha,mu)    
    implicit none
    real(kp) :: rclfi1_epsilon_three
    real(kp), intent(in) :: x,p,alpha,mu

    rclfi1_epsilon_three = rclfi_epsilon_three(x,p,alpha,mu)
    
  end function rclfi1_epsilon_three


  
!non vanishing x at which the potential is extremal
  function rclfi1_x_potmax(p,alpha,mu)
    implicit none
    real(kp) :: rclfi1_x_potmax
    real(kp), intent(in) :: p,alpha,mu

    real(kp), dimension(2) :: xdVzero

    if (.not.rclfi1_check_params(p,alpha,mu)) then
       stop 'rclfi1_x_potmax: no maximum!'
    endif

    xdVzero = rclfi_x_derivpotzero(p,alpha,mu)

    rclfi1_x_potmax = xdVzero(1)
    
  end function rclfi1_x_potmax



!returns the larger value of xpotmax in the -1 Lambert branch
  function rclfi1_max_xpotmax(p)
    implicit none
    real(kp) :: rclfi1_max_xpotmax
    real(kp), intent(in) :: p

    rclfi1_max_xpotmax = exp((8._kp-p)/(4._kp*(p-4._kp)))

  end function rclfi1_max_xpotmax



!returns the larger than 4 value of p such that the maximum of the
!potential remains larger than xpotmaxmin
  function rclfi1_numacc_pmax()
    implicit none
    real(kp) :: rclfi1_numacc_pmax

    rclfi1_numacc_pmax = (8._kp+16._kp*log(xpotmaxmin))/(1._kp+4._kp*log(xpotmaxmin))
    
  end function rclfi1_numacc_pmax
  
  

!return the value of alpha for which xpotmax=x in the branch where it
!can become arbitrarily small, i.e., for p<4 and alpha>alpha_one
  function rclfi1_numacc_alpha_xpotmax(x,p)
    implicit none
    real(kp) :: rclfi1_numacc_alpha_xpotmax
    real(kp), intent(in) :: x,p
    
    if (p.ge.4) then
       stop 'rclfi1_numacc_alpha_xpotmax: p>=4, incorrect branch!'
    endif

    if (x.gt.rclfi1_max_xpotmax(p)) then
       write(*,*)'rclfi1_numacc_alpha_xpotmax:'
       write(*,*)'x= max(xpotmax)= ',x,rclfi1_max_xpotmax(p)
       stop
    endif
    
    rclfi1_numacc_alpha_xpotmax = p*(p-4._kp)*x**(p-4._kp) &
         /(4._kp - p - 4._kp*(p-4._kp)*log(x))
    
    
  end function rclfi1_numacc_alpha_xpotmax
  

  function rclfi1_numacc_alphamax(p)
    implicit none
    real(kp) :: rclfi1_numacc_alphamax
    real(kp), intent(in) :: p
    
    rclfi1_numacc_alphamax = rclfi1_numacc_alpha_xpotmax(xpotmaxmin,p)
    
  end function rclfi1_numacc_alphamax
  


!returns x at the end of inflation defined as epsilon1=1
  function rclfi1_x_endinf(p,alpha,mu)
    implicit none
    real(kp), intent(in) :: p,alpha,mu
    real(kp) :: rclfi1_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rclfiData

    mini = epsilon(1._kp)
    maxi = rclfi1_x_potmax(p,alpha,mu)- epsilon(1._kp)

    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    
    rclfi1_x_endinf = zbrent(find_rclfi_x_endinf,mini,maxi,tolFind,rclfiData)
   
  end function rclfi1_x_endinf




  
!the closest possible to the top of the potential ensuring (dlnVodx)^2 > epsnumacc
  function rclfi1_numacc_xinimax(p,alpha,mu)
    implicit none
    real(kp) :: rclfi1_numacc_xinimax
    real(kp), intent(in) :: p,alpha,mu

    real(kp), parameter :: dx = 100._kp*epsnumacc
    real(kp) :: xpotmax, dlnVodx
    
    if (.not.rclfi1_check_params(p,alpha,mu)) then
       stop 'rclfi1_numacc_xinimax: no maximum!'
    endif

    xpotmax = rclfi1_x_potmax(p,alpha,mu)
    
    if (rclfi1_norm_potential(xpotmax-dx,p,alpha,mu).lt.0._kp) then
       write(*,*)'alpha= p= ',p,alpha
       stop 'rclfi1_numacc_xinimax: dx too large!'
    endif

    dlnVodx = abs(rclfi1_norm_deriv_potential(xpotmax-dx,p,alpha,mu) &
         /rclfi1_norm_potential(xpotmax,p,alpha,mu)/dx)
    
    rclfi1_numacc_xinimax = xpotmax - sqrt(epsnumacc)/dlnVodx

    
  end function rclfi1_numacc_xinimax



!maximal number of efolds computable at current numerical accuracy
  function rclfi1_numacc_efoldmax(p,alpha,mu)
    implicit none
    real(kp) :: rclfi1_numacc_efoldmax
    real(kp), intent(in) :: p,alpha,mu

    real(kp) :: xend,xinimax

    if (.not.rclfi1_check_params(p,alpha,mu)) then
       write(*,*)'alpha= p= ',p,alpha
       stop 'rclfi1_numacc_efoldmax: potential has no maximum!'
    endif

    xend = rclfi1_x_endinf(p,alpha,mu)

    xinimax = rclfi1_numacc_xinimax(p,alpha,mu)
    
    rclfi1_numacc_efoldmax = -rclfi1_efold_primitive(xend,p,alpha,mu) &
         + rclfi1_efold_primitive(xinimax,p,alpha,mu)
    
  end function rclfi1_numacc_efoldmax



!given alpha, p what is the minimal value of mu to get efold inflation
!above numerical accuracy limit
  function rclfi1_numacc_mumin(efold,p,alpha)
    implicit none
    real(kp) :: rclfi1_numacc_mumin
    real(kp), intent(in) :: efold,p,alpha

    real(kp), parameter :: tolFind=tolkp
    type(transfert) :: rclfiData

    real(kp), parameter :: mubig = 10000._kp
    real(kp), parameter :: musmall = 0.001_kp
    
    real(kp) :: mini,maxi

    mini = musmall
    maxi = mubig
    
    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = efold

    if (rclfi1_numacc_efoldmax(p,alpha,mubig).lt.efold) then
       stop 'rclfi1_numacc_mumin: mubig too small!'
    endif
    
    rclfi1_numacc_mumin = zbrent(find_rclfi1_numacc_mumin,mini,maxi,tolFind,rclfiData)

  end function rclfi1_numacc_mumin



  function find_rclfi1_numacc_mumin(mu,rclfiData)
    implicit none
    real(kp) :: find_rclfi1_numacc_mumin
    real(kp), intent(in) :: mu
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: alpha, p, Nwanted

    alpha = rclfiData%real1
    p = rclfiData%real2
    Nwanted = rclfiData%real3
        
    find_rclfi1_numacc_mumin = rclfi1_numacc_efoldmax(p,alpha,mu) - Nwanted
    
  end function find_rclfi1_numacc_mumin

  

!this is integral[V(phi)/V'(phi) dphi]
  function rclfi1_efold_primitive(x,p,alpha,mu)
    implicit none
    real(kp), intent(in) :: x,p,alpha,mu
    real(kp) :: rclfi1_efold_primitive

    type(transfert) :: rclfiData
!too long to integrate in QUADPREC
    real(kp), parameter :: tolInt = 10._kp*max(tolkp,toldp)
    integer, parameter :: neq = 1

    real(kp) :: xpotmax
    real(kp) :: xvar
    real(kp), dimension(neq) :: yvar


    if (x.gt.rclfi1_numacc_xinimax(p,alpha,mu)) then
       write(*,*)'rclfi1_efold_primitive: xVmax-x too small!'
       write(*,*)'x= alpha= p= mu= ',x,p,alpha,mu
       stop
    endif
    
    xpotmax = rclfi1_x_potmax(p,alpha,mu)

    if ((x.ge.xpotmax).or.(x.lt.0._kp)) then
       stop 'rclfi1_efold_primitive: x not in [0,xpotmax[!'
    endif

    xvar = 0._kp
    yvar(1) = 0._kp
    
    rclfiData%real1 = alpha
    rclfiData%real2 = p
    
    call easydverk(neq,find_rclfi_efold_primitive,xvar,yvar,x,tolInt,rclfiData)

    rclfi1_efold_primitive = yvar(1) * mu * mu
    
  end function rclfi1_efold_primitive


  
  
!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function rclfi1_x_trajectory(bfold,xend,p,alpha,mu)
    implicit none
    real(kp), intent(in) :: bfold, alpha, p, mu, xend
    real(kp) :: rclfi1_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rclfiData

  
    mini = xend
    maxi = rclfi1_numacc_xinimax(p,alpha,mu)
  
    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    rclfiData%real4 = -bfold + rclfi1_efold_primitive(xend,p,alpha,mu)

    rclfi1_x_trajectory = zbrent(find_rclfi1_x_trajectory,mini,maxi,tolFind,rclfiData)
    
  end function rclfi1_x_trajectory


  
  function find_rclfi1_x_trajectory(x,rclfiData)
    implicit none
    real(kp) :: find_rclfi1_x_trajectory
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: alpha, p, mu, NplusNuend

    alpha = rclfiData%real1
    p = rclfiData%real2
    mu = rclfiData%real3    
    NplusNuend = rclfiData%real4

    find_rclfi1_x_trajectory = rclfi1_efold_primitive(x,p,alpha,mu) - NplusNuend

  end function find_rclfi1_x_trajectory

  
end module rclfi1sr
