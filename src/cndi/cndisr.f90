!slow-roll functions for the non canonical Kahler inflation potential
!
!V(phi) = M**4 / ( 1 + beta cos(alpha x) )**2
!
!x = (phi-phi0)/Mp

module cndisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public cndi_norm_potential, cndi_epsilon_one, cndi_epsilon_two, cndi_epsilon_three
  public cndi_efold_primitive, cndi_x_trajectory
  public cndi_norm_deriv_potential, cndi_norm_deriv_second_potential
  public cndi_xinimax, cndi_xendmax

 
contains
!returns V/M**4 as function of x
  function cndi_norm_potential(x,alpha,beta)
    implicit none
    real(kp) :: cndi_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    cndi_norm_potential = 1._kp/(1._kp+beta*cos(alpha*x))**2

  end function cndi_norm_potential



!returns the first derivative of the potential with respect to x
  function cndi_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: cndi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta

   cndi_norm_deriv_potential = (2._kp*alpha*beta*sin(alpha*x))/ &
        (1._kp+beta*cos(alpha*x))**3

  end function cndi_norm_deriv_potential



!returns the second derivative of the potential with respect to x
  function cndi_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: cndi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta

    cndi_norm_deriv_second_potential = (2._kp*alpha**2*beta*(cos(alpha*x)- &
         beta*(-2._kp+cos(2._kp*alpha*x))))/ & 
         (1._kp+beta*cos(alpha*x))**4

  end function cndi_norm_deriv_second_potential



!epsilon_one(x)
  function cndi_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: cndi_epsilon_one
    real(kp), intent(in) :: x,alpha,beta
    
    cndi_epsilon_one = (2._kp*alpha**2*beta**2*sin(alpha*x)**2)/ &
         (1._kp+beta*cos(alpha*x))**2
    
  end function cndi_epsilon_one


!epsilon_two(x)
  function cndi_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: cndi_epsilon_two
    real(kp), intent(in) :: x,alpha,beta
    
    cndi_epsilon_two = -((4._kp*alpha**2*beta*(beta+cos(alpha*x)))/ &
         (1._kp+beta*cos(alpha*x))**2)
    
  end function cndi_epsilon_two


!epsilon_three(x)
  function cndi_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: cndi_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
    
    cndi_epsilon_three = -((2._kp*alpha**2*beta*(-1._kp+2._kp*beta**2+ &
         beta*cos(alpha*x))*sin(alpha*x)**2)/((beta+cos(alpha*x))* &
         (1._kp+beta*cos(alpha*x))**2))
    
  end function cndi_epsilon_three



!returns the minimum value of alpha, given beta, such that eps1=1 has a solution for beta<1
  function cndi_alphamin(beta)
    implicit none
    real(kp) , intent(in) :: beta
    real(kp) :: cndi_alphamin

    if (beta.gt.1._kp) then
       cndi_alphamin=0._kp
    else
       cndi_alphamin=sqrt((1._kp-beta**2)/(2._kp*beta**2))
    endif

  end function cndi_alphamin



  function cndi_x_epsoneunity(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp), dimension(2) :: cndi_x_epsoneunity

    if (alpha.lt.cndi_alphamin(beta)) then
       print*,'alpha=',alpha,'beta=',beta,'alphamin=',cndi_alphamin(beta)
       stop 'cndi_x_epsoneunity: not root for epsilon1 = 1' 
    endif

    cndi_x_epsoneunity(1) = acos((alpha*sqrt(beta**2*(2._kp+4._kp*alpha**2)-2._kp)-1._kp)/ &
         (beta+2._kp*alpha**2*beta))/alpha

    cndi_x_epsoneunity(2) = acos(-((alpha*sqrt(beta**2*(2._kp+4._kp*alpha**2)-2._kp)+1._kp))/ &
         (beta+2._kp*alpha**2*beta))/alpha

  end function cndi_x_epsoneunity




!returns the maximum authorized value for x (in order for epsilon1<1)
  function cndi_xinimax(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: cndi_xinimax

    real(kp), dimension(2) :: xepsone
    
    if (beta .gt. 1._kp .or. (beta .lt. 1._kp .and. alpha.gt.cndi_alphamin(beta) ) ) then
      xepsone = cndi_x_epsoneunity(alpha,beta)
      cndi_xinimax  = xepsone(1)
    else
      cndi_xinimax  = acos(-1._kp)/alpha*(1._kp-epsilon(1._kp))
    endif
   
  end function cndi_xinimax

!returns the maximum authorized value for xend in order for efold to be realized from xinimax
  function cndi_xendmax(efold,alpha,beta)
    implicit none
    real(kp), intent(in) :: efold,alpha,beta
    real(kp) :: cndi_xendmax,xiniMax,xendMin

    real(kp) :: efoldNumAccMax

    xiniMax = cndi_xinimax(alpha,beta)
    xendMin = epsilon(1._kp)

    efoldNumAccMax = -cndi_efold_primitive(xendMin,alpha,beta) &
         + cndi_efold_primitive(xiniMax,alpha,beta)

    if (efold.gt.efoldNumAccMax) then
       write(*,*)'cndi_xendmax: not enough efolds computable at current accuracy!'
       write(*,*)'efold requested=   efold maxi computable= ',efold,efoldNumAccMax
       write(*,*)'xend mini computable= ',xendMin
       stop
    endif
   
    cndi_xendmax = cndi_x_trajectory(efold,xiniMax,alpha,beta)

   
  end function cndi_xendmax


!this is integral(V(phi)/V'(phi) dphi)
  function cndi_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: cndi_efold_primitive

    if (alpha*beta.eq.0._kp) stop 'cndi_efold_primitive: alpha*beta=0 !'

    cndi_efold_primitive = (log(sin(alpha*x))+log(tan(0.5_kp*alpha*x))/beta)/(2._kp*alpha**2)

  end function cndi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function cndi_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold,alpha,beta,xend
    real(kp) :: cndi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: cndiData

    cndiData%real1 = alpha
    cndiData%real2 = beta
    cndiData%real3 = -bfold + cndi_efold_primitive(xend,alpha,beta)


    if (bfold .lt. 0._kp) then
      maxi = cndi_xinimax(alpha,beta)
      mini = xend
    else
      maxi = xend
!to avoid < numaccuracy errors
!      mini = sqrt(epsilon(1._kp)/2._kp)*(1._kp+beta)/(alpha**2*beta) 
      mini = epsilon(1._kp)
    endif
    
    cndi_x_trajectory = zbrent(find_cndi_x_trajectory,mini,maxi,tolFind,cndiData)
       
  end function cndi_x_trajectory

  function find_cndi_x_trajectory(x,cndiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: cndiData
    real(kp) :: find_cndi_x_trajectory
    real(kp) :: alpha,beta,NplusNuend

    alpha= cndiData%real1
    beta = cndiData%real2
    NplusNuend = cndiData%real3

    find_cndi_x_trajectory = cndi_efold_primitive(x,alpha,beta) - NplusNuend
   
  end function find_cndi_x_trajectory


  
end module cndisr
