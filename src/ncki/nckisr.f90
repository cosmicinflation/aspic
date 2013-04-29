!slow-roll functions for the non canonical Kahler inflation potential
!
!V(phi) = M^4 [ 1 + alpha ln(x) + beta (x)^2 ]
!
!x = phi/Mp

module nckisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : polylog,lambert
  implicit none

  private

  real(kp), parameter :: NckiMaxiInf = 1000._kp

  public NckiMaxiInf
  public ncki_norm_potential, ncki_epsilon_one, ncki_epsilon_two, ncki_epsilon_three
  public ncki_x_endinf, ncki_efold_primitive, ncki_x_trajectory
  public ncki_norm_deriv_potential, ncki_norm_deriv_second_potential
  public ncki_x_potmax, ncki_x_inflection, ncki_x_dderivpotzero


contains
  !returns V/M**4 as function of x
  function ncki_norm_potential(x,alpha,beta)
    implicit none
    real(kp) :: ncki_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ncki_norm_potential = 1._kp+alpha*log(x)+beta*x**2

  end function ncki_norm_potential



  !returns the first derivative of the potential with respect to x
  function ncki_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ncki_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta

    ncki_norm_deriv_potential = alpha/x+2._kp*beta*x

  end function ncki_norm_deriv_potential



  !returns the second derivative of the potential with respect to x
  function ncki_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ncki_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta

    ncki_norm_deriv_second_potential = 2._kp*beta-alpha/(x**2)

  end function ncki_norm_deriv_second_potential



  !epsilon_one(x)
  function ncki_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ncki_epsilon_one
    real(kp), intent(in) :: x,alpha,beta

    ncki_epsilon_one = (alpha+2._kp*beta*x**2)**2/(2._kp*(x+beta*x**3+alpha*x*log(x))**2)

  end function ncki_epsilon_one


  !epsilon_two(x)
  function ncki_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ncki_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    ncki_epsilon_two = (2._kp*(alpha*(1._kp+alpha)+(-2._kp+5._kp*alpha)*beta*x**2 &
         +2._kp*beta**2*x**4+alpha*(alpha-2._kp*beta*x**2)*log(x)))/ &
         (x+beta*x**3+alpha*x*log(x))**2

  end function ncki_epsilon_two


  !epsilon_three(x)
  function ncki_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ncki_epsilon_three
    real(kp), intent(in) :: x,alpha,beta

    ncki_epsilon_three = (1._kp/(x**2))*((2._kp*(alpha+2._kp*beta*x**2)**2)/(1._kp+beta*x**2 &
         +alpha*log(x))**2+(alpha-2._kp*beta*x**2)/(1._kp+beta*x**2+alpha*log(x)) &
         +(alpha**2+8._kp*alpha*beta*x**2-4._kp*beta**2*x**4)/(alpha*(1._kp+alpha) & 
         +(-2._kp+5._kp*alpha)*beta*x**2+2._kp*beta**2*x**4+alpha*(alpha-2._kp*beta*x**2)*log(x)))

  end function ncki_epsilon_three



!field value at which the potential is maximal when beta<0
  function ncki_x_potmax(alpha,beta)
    implicit none
    real(kp) , intent(in) :: alpha,beta
    real(kp) :: ncki_x_potmax

    if (beta.gt.0._kp) stop 'ncki_x_potmax: no maximum for beta>0!'

    ncki_x_potmax = ncki_x_dderivpotzero(alpha,beta)

  end function ncki_x_potmax


!field value at which the potential is maximal when beta>0
  function ncki_x_inflection(alpha,beta)
    implicit none
    real(kp) , intent(in) :: alpha,beta
    real(kp) :: ncki_x_inflection

    if (beta.lt.0._kp) stop 'ncki_x_inflection: no inflection for beta<0!'

    ncki_x_inflection = ncki_x_dderivpotzero(alpha,beta)

  end function ncki_x_inflection


  function ncki_x_dderivpotzero(alpha,beta)
    implicit none
    real(kp) , intent(in) :: alpha,beta
    real(kp) :: ncki_x_dderivpotzero

    ncki_x_dderivpotzero = sqrt(abs(alpha/2._kp/beta))

  end function ncki_x_dderivpotzero


  !returns x at the end of inflation defined as epsilon1=1
  function ncki_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: ncki_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: nckiData

!mini =
!sqrt(alpha/(2._kp*beta)*lambert(2._kp*beta/alpha*exp(-2._kp/alpha),0))*(1._kp+epsilon(1._kp))
!minimum value below which the potential is negative
    mini=epsilon(1._kp)
!position of the maximum of the potential if beta<0, and of the inflexion point if beta>0
    maxi = ncki_x_dderivpotzero(alpha,beta)*(1._kp-epsilon(1._kp))

    nckiData%real1 = alpha
    nckiData%real2 = beta

    ncki_x_endinf = zbrent(find_ncki_x_endinf,mini,maxi,tolFind,nckiData)

  end function ncki_x_endinf

  function find_ncki_x_endinf(x,nckiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: nckiData
    real(kp) :: find_ncki_x_endinf
    real(kp) :: alpha,beta

    alpha = nckiData%real1
    beta = nckiData%real2

    find_ncki_x_endinf = ncki_epsilon_one(x,alpha,beta)-1._kp

  end function find_ncki_x_endinf


  !this is integral(V(phi)/V'(phi) dphi)
  function ncki_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ncki_efold_primitive

    if (alpha.eq.0._kp) stop 'ncki_efold_primitive: alpha=0 !'

    ncki_efold_primitive = (1._kp-alpha/2._kp+alpha*log(x)) &
         *log(alpha+2._kp*beta*x**2)/(4._kp*beta) &
         +x**2/4._kp-alpha/(4._kp*beta)*log(alpha)*log(x)+alpha/(8._kp*beta)* &
         real(polylog(complex(-2._kp*beta/alpha*x**2,0._kp),complex(2._kp,0._kp)),kp)


  end function ncki_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ncki_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold,alpha,beta,xend
    real(kp) :: ncki_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: nckiData

    if (beta.lt.0._kp) then
!position of the maximum of the potential if beta<0
       maxi =ncki_x_dderivpotzero(alpha,beta)*(1._kp-epsilon(1._kp))
    else
!several times the position of the inflexion point if beta>0
       maxi = NckiMaxiInf * ncki_x_dderivpotzero(alpha,beta)
    endif

    mini=epsilon(1._kp)

    nckiData%real1 = alpha
    nckiData%real2 = beta
    nckiData%real3 = -bfold + ncki_efold_primitive(xend,alpha,beta)

    ncki_x_trajectory = zbrent(find_ncki_x_trajectory,mini,maxi,tolFind,nckiData)

  end function ncki_x_trajectory

  function find_ncki_x_trajectory(x,nckiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: nckiData
    real(kp) :: find_ncki_x_trajectory
    real(kp) :: alpha,beta,NplusNuend

    alpha= nckiData%real1
    beta = nckiData%real2
    NplusNuend = nckiData%real3

    find_ncki_x_trajectory = ncki_efold_primitive(x,alpha,beta) - NplusNuend

  end function find_ncki_x_trajectory



end module nckisr
