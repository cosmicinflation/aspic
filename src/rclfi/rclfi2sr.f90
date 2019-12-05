!Radiatively corrected large field inflation in the intermediate field
!region; after the first maximum of the potential and before the
!second minimum (required to be <=0)
!
!V(phi) = M^4 [x^p + alpha x^4 ln(x) ]
!
!x = phi/mu
!
!with:
!  p>4 and alpha< -e(p-4) [< -p(p-4)/4 exp(2-p/4)]
!OR
!  p<4 and alpha<0
!OR
!  p<4 and alpha> -e(p-4) [>-p(p-4)/4 exp(2-p/4)]
!

module rclfi2sr
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

  public rclfi2_check_params, rclfi2_norm_deriv_second_potential
  public rclfi2_norm_potential, rclfi2_norm_deriv_potential
  public rclfi2_epsilon_one, rclfi2_epsilon_two, rclfi2_epsilon_three
  public rclfi2_x_endinf, rclfi2_efold_primitive, rclfi2_x_trajectory
  public rclfi2_x_potmax, rclfi2_numacc_xinimin, rclfi2_numacc_efoldmax
  public rclfi2_numacc_mumin, rclfi2_numacc_pmax, rclfi2_numacc_alphamax

!numerical accuracy limitation
  real(kp), parameter :: epsnumacc = 100._kp*tolkp

!we do not consider inflation occuring at xpotmax < x< xpotzero when
!xpotzero-xpotmax < deltaxmin
  real(kp), parameter :: deltaxmin = 0.1_kp
  
contains


  
  function rclfi2_check_params(alpha,p,mu)
    logical :: rclfi2_check_params
    real(kp), intent(in) :: alpha,p,mu

    real(kp) :: alphaVo

    alphaVo = rclfi_alpha_zero(p)
    
    rclfi2_check_params = ((p.gt.4._kp).and.(alpha.le.alphaVo)) &
         .or. ((p.lt.4._kp).and.(alpha.lt.0._kp)) &
         .or. ((p.lt.4._kp).and.(alpha.ge.alphaVo))
    
    rclfi2_check_params = rclfi2_check_params.and.rclfi_check_params(alpha,p,mu)
    
  end function rclfi2_check_params
  

  
  function rclfi2_norm_potential(x,alpha,p,mu)
    implicit none
    real(kp) :: rclfi2_norm_potential
    real(kp), intent(in) :: x,alpha,p,mu

    rclfi2_norm_potential = rclfi_norm_potential(x,alpha,p,mu)

  end function rclfi2_norm_potential




  function rclfi2_norm_deriv_potential(x,alpha,p,mu)
    implicit none
    real(kp) :: rclfi2_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,p,mu

   rclfi2_norm_deriv_potential = rclfi_norm_deriv_potential(x,alpha,p,mu)

  end function rclfi2_norm_deriv_potential




  function rclfi2_norm_deriv_second_potential(x,alpha,p,mu)
    implicit none
    real(kp) :: rclfi2_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,p,mu

    rclfi2_norm_deriv_second_potential = rclfi_norm_deriv_second_potential(x,alpha,p,mu)

  end function rclfi2_norm_deriv_second_potential




  function rclfi2_epsilon_one(x,alpha,p,mu)    
    implicit none
    real(kp) :: rclfi2_epsilon_one
    real(kp), intent(in) :: x,alpha,p,mu

    rclfi2_epsilon_one = rclfi_epsilon_one(x,alpha,p,mu)
    
  end function rclfi2_epsilon_one



  function rclfi2_epsilon_two(x,alpha,p,mu)    
    implicit none
    real(kp) :: rclfi2_epsilon_two
    real(kp), intent(in) :: x,alpha,p,mu

    rclfi2_epsilon_two = rclfi_epsilon_two(x,alpha,p,mu)
    
  end function rclfi2_epsilon_two



  function rclfi2_epsilon_three(x,alpha,p,mu)    
    implicit none
    real(kp) :: rclfi2_epsilon_three
    real(kp), intent(in) :: x,alpha,p,mu

    rclfi2_epsilon_three = rclfi_epsilon_three(x,alpha,p,mu)
    
  end function rclfi2_epsilon_three


  
!non vanishing x at which the potential is extremal
  function rclfi2_x_potmax(alpha,p,mu)
    implicit none
    real(kp) :: rclfi2_x_potmax
    real(kp), intent(in) :: alpha,p,mu

    real(kp), dimension(2) :: xdVzero

    if (.not.rclfi2_check_params(alpha,p,mu)) then
       stop 'rclfi2_x_potmax: no maximum!'
    endif

    xdVzero = rclfi_x_derivpotzero(alpha,p,mu)

    rclfi2_x_potmax = xdVzero(1)
    
  end function rclfi2_x_potmax


!non vanishing x at which the potential becomes negative or null
  function rclfi2_x_potzero(alpha,p,mu)
    implicit none
    real(kp) :: rclfi2_x_potzero
    real(kp), intent(in) :: alpha,p,mu

    real(kp), dimension(2) :: xVzero

    if (.not.rclfi2_check_params(alpha,p,mu)) then
       stop 'rclfi2_x_potzero: no zeros!'
    endif
    
    xVzero = rclfi_x_potzero(alpha,p,mu)

    rclfi2_x_potzero = xVzero(1)
    
  end function rclfi2_x_potzero

  

!returns the larger value of xpotmax in the -1 Lambert branch
  function rclfi2_max_xpotmax(p)
    implicit none
    real(kp) :: rclfi2_max_xpotmax
    real(kp), intent(in) :: p

    rclfi2_max_xpotmax = exp((8._kp-p)/(4._kp*(p-4._kp)))

  end function rclfi2_max_xpotmax


  

!returns x at the end of inflation defined as epsilon1=1
  function rclfi2_x_endinf(alpha,p,mu)
    implicit none
    real(kp), intent(in) :: alpha,p,mu
    real(kp) :: rclfi2_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rclfiData

    mini = rclfi2_x_potmax(alpha,p,mu)+ epsilon(1._kp)
    maxi = rclfi2_x_potzero(alpha,p,mu)- epsilon(1._kp)

    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    
    rclfi2_x_endinf = zbrent(find_rclfi_x_endinf,mini,maxi,tolFind,rclfiData)
   
  end function rclfi2_x_endinf


  
!see next function, that's the pmax numerically achievable  
  function rclfi2_numacc_pmax()
    implicit none
    real(kp) :: rclfi2_numacc_pmax
    
    rclfi2_numacc_pmax = rclfi2_numacc_pmax_deltax(deltaxmin)
    
  end function rclfi2_numacc_pmax
  

  

!return the maximal value of p at alphamin=-e(p-4) for which
!xpotzero-xpotmax=x in the branch where it can become arbitrarily
!small, i.e., for p<4 
  function rclfi2_numacc_pmax_deltax(deltax)
    implicit none
    real(kp) :: rclfi2_numacc_pmax_deltax
    real(kp), intent(in) :: deltax
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rclfiData

    real(kp), parameter :: junk = -1_kp
              
       
    mini = 100*tolkp
    maxi = 4._kp - 100*tolkp

    rclfiData%real1=deltax
    
    rclfi2_numacc_pmax_deltax = zbrent(find_rclfi2_numacc_pmax_deltax,mini,maxi,tolFind,rclfiData)
        
  end function rclfi2_numacc_pmax_deltax

  

  function find_rclfi2_numacc_pmax_deltax(p,rclfiData)
    implicit none
    real(kp) :: find_rclfi2_numacc_pmax_deltax
    real(kp), intent(in) :: p
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: deltax,alpha
    real(kp), parameter :: junk = -1._kp

    deltax = rclfiData%real1

    alpha = -exp(1._kp)*(p-4._kp)*(1._kp + tolkp)

    find_rclfi2_numacc_pmax_deltax = rclfi2_x_potzero(alpha,p,junk) &
         - rclfi2_x_potmax(alpha,p,junk) - deltax

  end function find_rclfi2_numacc_pmax_deltax



  
!see next one  
  function rclfi2_numacc_alphamax(p)
    implicit none
    real(kp) :: rclfi2_numacc_alphamax
    real(kp), intent(in) :: p

    rclfi2_numacc_alphamax = rclfi2_numacc_alphamax_deltax(deltaxmin,p)
        
  end function rclfi2_numacc_alphamax

  
  
!return the maximal value of alpha at p<pmax<4 for which
!xpotzero-xpotmax=x in the branch where it can become arbitrarily
!small
  function rclfi2_numacc_alphamax_deltax(deltax,p)
    implicit none
    real(kp) :: rclfi2_numacc_alphamax_deltax
    real(kp), intent(in) :: deltax,p
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rclfiData

    real(kp), parameter :: junk = -1_kp
              
       
    mini = -exp(1._kp)*(p-4._kp) + tolkp
    maxi = 1._kp/epsilon(1._kp)

    rclfiData%real1=deltax
    rclfiData%real2=p
    
    rclfi2_numacc_alphamax_deltax &
         = zbrent(find_rclfi2_numacc_alphamax_deltax,mini,maxi,tolFind,rclfiData)
        
  end function rclfi2_numacc_alphamax_deltax

  

  function find_rclfi2_numacc_alphamax_deltax(alpha,rclfiData)
    implicit none
    real(kp) :: find_rclfi2_numacc_alphamax_deltax
    real(kp), intent(in) :: alpha
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: deltax,p
    real(kp), parameter :: junk = -1._kp

    deltax = rclfiData%real1
    p = rclfiData%real2

    find_rclfi2_numacc_alphamax_deltax = rclfi2_x_potzero(alpha,p,junk) &
         - rclfi2_x_potmax(alpha,p,junk) - deltax

  end function find_rclfi2_numacc_alphamax_deltax

  
  
!the closest possible to the top of the potential ensuring (dlnVodx)^2 > epsnumacc
  function rclfi2_numacc_xinimin(alpha,p,mu)
    implicit none
    real(kp) :: rclfi2_numacc_xinimin
    real(kp), intent(in) :: alpha,p,mu

    real(kp), parameter :: dx = 10._kp*epsnumacc
    real(kp) :: xpotmax, dlnVodx
    
    if (.not.rclfi2_check_params(alpha,p,mu)) then
       stop 'rclfi2_numacc_xinimax: no maximum!'
    endif

    xpotmax = rclfi2_x_potmax(alpha,p,mu)
    
    if (rclfi2_norm_potential(xpotmax+dx,alpha,p,mu).lt.0._kp) then
       write(*,*)'alpha= p= ',alpha,p
       stop 'rclfi2_numacc_xinimin: dx too large!'
    endif

    dlnVodx = abs(rclfi2_norm_deriv_potential(xpotmax+dx,alpha,p,mu) &
         /rclfi2_norm_potential(xpotmax,alpha,p,mu)/dx)
    
    rclfi2_numacc_xinimin = xpotmax + sqrt(epsnumacc)/dlnVodx

    
  end function rclfi2_numacc_xinimin



!maximal number of efolds computable at current numerical accuracy
  function rclfi2_numacc_efoldmax(alpha,p,mu)
    implicit none
    real(kp) :: rclfi2_numacc_efoldmax
    real(kp), intent(in) :: alpha,p,mu

    real(kp) :: xend,xinimin

    if (.not.rclfi2_check_params(alpha,p,mu)) then
       write(*,*)'alpha= p= ',alpha,p
       stop 'rclfi2_numacc_efoldmax: potential has no maximum!'
    endif

    xend = rclfi2_x_endinf(alpha,p,mu)

    xinimin = rclfi2_numacc_xinimin(alpha,p,mu)
    
    rclfi2_numacc_efoldmax = -rclfi2_efold_primitive(xend,alpha,p,mu) &
         + rclfi2_efold_primitive(xinimin,alpha,p,mu)
    
  end function rclfi2_numacc_efoldmax



!given alpha, p what is the minimal value of mu to get efold inflation
!above numerical accuracy limit
  function rclfi2_numacc_mumin(efold,alpha,p)
    implicit none
    real(kp) :: rclfi2_numacc_mumin
    real(kp), intent(in) :: efold,alpha,p

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

    if (rclfi2_numacc_efoldmax(alpha,p,mubig).lt.efold) then
       write(*,*)'max efold= ',rclfi2_numacc_efoldmax(alpha,p,mubig)
       stop 'rclfi2_numacc_mumin: mubig too small!'
    endif
    
    rclfi2_numacc_mumin = zbrent(find_rclfi2_numacc_mumin,mini,maxi,tolFind,rclfiData)

  end function rclfi2_numacc_mumin



  function find_rclfi2_numacc_mumin(mu,rclfiData)
    implicit none
    real(kp) :: find_rclfi2_numacc_mumin
    real(kp), intent(in) :: mu
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: alpha, p, Nwanted

    alpha = rclfiData%real1
    p = rclfiData%real2
    Nwanted = rclfiData%real3
        
    find_rclfi2_numacc_mumin = rclfi2_numacc_efoldmax(alpha,p,mu) - Nwanted
    
  end function find_rclfi2_numacc_mumin

  

!this is integral[V(phi)/V'(phi) dphi]
  function rclfi2_efold_primitive(x,alpha,p,mu)
    implicit none
    real(kp), intent(in) :: x,alpha,p,mu
    real(kp) :: rclfi2_efold_primitive

    type(transfert) :: rclfiData
!too long to integrate in QUADPREC
    real(kp), parameter :: tolInt = 10._kp*max(tolkp,toldp)
    integer, parameter :: neq = 1

    real(kp) :: xpotmax,xpotzero
    real(kp) :: xvar
    real(kp), dimension(neq) :: yvar


    if (x.lt.rclfi2_numacc_xinimin(alpha,p,mu)) then
       write(*,*)'rclfi2_efold_primitive: xVmax+x too small!'
       write(*,*)'x= alpha= p= mu= ',x,alpha,p,mu
       stop
    endif
    
    xpotmax = rclfi2_x_potmax(alpha,p,mu)
    xpotzero = rclfi2_x_potzero(alpha,p,mu)
    
    if ((x.le.xpotmax).or.(x.ge.xpotzero)) then
       stop 'rclfi2_efold_primitive: x not in ]xpotmax,xpotzero[!'
    endif

    xvar = xpotzero
    yvar(1) = 0._kp
    
    rclfiData%real1 = alpha
    rclfiData%real2 = p
    
    call easydverk(neq,find_rclfi_efold_primitive,xvar,yvar,x,tolInt,rclfiData)

    rclfi2_efold_primitive = yvar(1) * mu * mu
    
  end function rclfi2_efold_primitive


  
  
!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function rclfi2_x_trajectory(bfold,xend,alpha,p,mu)
    implicit none
    real(kp), intent(in) :: bfold, alpha, p, mu, xend
    real(kp) :: rclfi2_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rclfiData


    mini = rclfi2_numacc_xinimin(alpha,p,mu)
    maxi = xend

  
    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    rclfiData%real4 = -bfold + rclfi2_efold_primitive(xend,alpha,p,mu)

    rclfi2_x_trajectory = zbrent(find_rclfi2_x_trajectory,mini,maxi,tolFind,rclfiData)
    
  end function rclfi2_x_trajectory


  
  function find_rclfi2_x_trajectory(x,rclfiData)
    implicit none
    real(kp) :: find_rclfi2_x_trajectory
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: alpha, p, mu, NplusNuend

    alpha = rclfiData%real1
    p = rclfiData%real2
    mu = rclfiData%real3    
    NplusNuend = rclfiData%real4

    find_rclfi2_x_trajectory = rclfi2_efold_primitive(x,alpha,p,mu) - NplusNuend

  end function find_rclfi2_x_trajectory

  
end module rclfi2sr
