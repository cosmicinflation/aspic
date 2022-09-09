!slow-roll functions for the orientifold inflation potential
!
!V(phi) = M**4 / [ 1 + alpha log(phi/Mp) + (phi/phi0)**(4+2n) ]
!
!x = phi/Mp

module nclisr
  use infprec, only : kp, toldp, tolkp, transfert
  use inftools, only : zbrent, easydverk
  use specialinf, only : lambert
  implicit none

  private


!maximal value of x numerically usable
  real(kp), parameter :: xNumAccMax = 1000._kp



  public ncli_norm_potential, ncli_norm_deriv_potential, ncli_norm_deriv_second_potential
  public ncli_epsilon_one, ncli_epsilon_two,ncli_epsilon_three
  public ncli_efold_primitive, ncli_x_trajectory, ncli_x_endinf, ncli_check_params
  public ncli_x_epstwozero, ncli_x_epsoneunity, ncli_epsilon_one_max, ncli_epsilon_one_min
  public ncli_xinimin, ncli_xinimax, ncli_x_potzero, ncli_x_inflection, ncli_phizeromin
  public xNumAccMax



contains

  function ncli_check_params(efoldMin,epsMin,alpha,phi0,n)
    implicit none
    logical :: ncli_check_params
    real(kp), intent(in) :: efoldMin,epsMin
    real(kp), intent(in) :: alpha,phi0,n

    real(kp) :: eps1Min, eps1Max, xiniMax
    real(kp) :: efoldMax, xend

    ncli_check_params = .true.

    if (phi0.le.ncli_phizeromin(alpha,n)) then
       ncli_check_params = .false.
       return
    endif

    eps1min = ncli_epsilon_one_min(alpha,phi0,n)
    eps1max = ncli_epsilon_one_max(alpha,phi0,n)

    if (eps1min.ge.epsMin) then
       ncli_check_params = .false.
       return
    endif

    xend = ncli_x_endinf(alpha,phi0,n)
    xiniMax = ncli_xinimax(alpha,phi0,n)

    efoldMax = -ncli_efold_primitive(xend,alpha,phi0,n) &
         + ncli_efold_primitive(xiniMax,alpha,phi0,n)

    if (efoldMax.lt.efoldMin) then
       ncli_check_params = .false.
    endif

    if (phi0.lt.ncli_phizeromin(alpha,n)) then
      ncli_check_params = .false.
    endif


  end function ncli_check_params


!returns V/M**4
  function ncli_norm_potential(x,alpha,phi0,n)
    implicit none
    real(kp) :: ncli_norm_potential
    real(kp), intent(in) :: x,alpha,phi0,n

    ncli_norm_potential = 1._kp + 0.5_kp*alpha*log(x*x) + (x*x/phi0/phi0)**(2._kp+n)

  end function ncli_norm_potential


  !returns the first derivative of the potential with respect to x=phi/phi0, divided by M**4
  function ncli_norm_deriv_potential(x,alpha,phi0,n)
    implicit none
    real(kp) :: ncli_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,phi0,n

    ncli_norm_deriv_potential = alpha/x+(2._kp*(2._kp+n)*x**3*(x*x/phi0/phi0)**n)/phi0**4

  end function ncli_norm_deriv_potential



  !returns the second derivative of the potential with respect to x, divided by M**4
  function ncli_norm_deriv_second_potential(x,alpha,phi0,n)
    implicit none
    real(kp) :: ncli_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,phi0,n

    ncli_norm_deriv_second_potential = -(alpha/x**2)+(2._kp*(2._kp+n)*(3._kp+2._kp*n)*x**2* &
         (x*x/phi0/phi0)**n)/phi0**4

  end function ncli_norm_deriv_second_potential


  function ncli_epsilon_one(x,alpha,phi0,n)
    implicit none
    real(kp) :: ncli_epsilon_one
    real(kp), intent(in) :: x,alpha,phi0,n

    ncli_epsilon_one =(alpha/x+(2._kp*(2._kp+n)*x**3*(x/phi0)**(2._kp*n))/phi0**4)**2/ &
         (2._kp*(1._kp+(x/phi0)**(2._kp*(2._kp+n))+alpha*log(x))**2)

  end function ncli_epsilon_one



  function ncli_epsilon_two(x,alpha,phi0,n)
    implicit none
    real(kp) :: ncli_epsilon_two
    real(kp), intent(in) :: x,alpha,phi0,n

    ncli_epsilon_two = 2._kp*((alpha/x+(2._kp*(2._kp+n)*x**3*(x/phi0)**(2._kp*n))/phi0**4)**2/ &
         (1._kp+(x/phi0)**(2._kp*(2._kp+n))+alpha*log(x))**2-(-(alpha/x**2)+ &
         (2._kp*(2._kp+n)*(3._kp+2._kp*n)*x**2*(x/phi0)**(2._kp*n))/phi0**4)/ &
         (1._kp+(x/phi0)**(2._kp*(2._kp+n))+alpha*log(x)))

  end function ncli_epsilon_two


  function ncli_epsilon_three(x,alpha,phi0,n)
    implicit none
    real(kp) :: ncli_epsilon_three
    real(kp), intent(in) :: x,alpha,phi0,n

    ncli_epsilon_three = ((x**5*(x/phi0)**(2._kp*n)+x*phi0**4*(1._kp+alpha*log(x)))**2* &
         ((2._kp*(alpha/x+(2._kp*(2._kp+n)*x**3*(x/phi0)**(2._kp*n))/ &
         phi0**4)**4)/(1._kp+(x/phi0)**(2._kp*(2._kp+n))+alpha*log(x))**4- &
         (3._kp*(-(alpha/x**2)+(2._kp*(2._kp+n)*(3._kp+2._kp*n)*x**2* &
         (x/phi0)**(2.*n))/phi0**4)*(alpha/x+(2._kp*(2._kp+n)*x**3* &
         (x/phi0)**(2._kp*n))/phi0**4)**2)/(1._kp+(x/phi0)**(2._kp*(2._kp+n))+ &
         alpha*log(x))**3+(((2._kp*alpha)/x**3+(4._kp*(1._kp+n)*(2._kp+n)* &
         (3._kp+2._kp*n)*x*(x/phi0)**(2._kp*n))/phi0**4)*(alpha/x+(2._kp* &
         (2._kp+n)*x**3*(x/phi0)**(2._kp*n))/phi0**4))/(1._kp+(x/phi0)**( &
         2._kp*(2._kp+n))+alpha*log(x))**2))/(2._kp*(2._kp+n)*x**8*(x/phi0)**( &
         4._kp*n)+alpha*phi0**8*(1._kp+alpha+alpha*log(x))-x**4*(x/phi0)**( &
         2._kp*n)*phi0**4*(12._kp+2._kp*n*(7._kp+2._kp*n-2._kp*alpha)-9._kp* &
         alpha+2._kp*(2._kp+n)*(3._kp+2._kp*n)*alpha*log(x)))

  end function ncli_epsilon_three



!returns the point where the potential vanishes
  function ncli_x_potzero(alpha,phi0,n)
    implicit none
    real(kp) :: ncli_x_potzero
    real(kp), intent(in) :: alpha,phi0,n

    ncli_x_potzero = phi0*(alpha/(4._kp+2._kp*n)*lambert((4._kp+2._kp*n)/alpha* &
         exp(-(4._kp+2._kp*n)*(log(phi0)+1._kp/alpha)),0))**(1._kp/(4._kp+2._kp*n))

  end function ncli_x_potzero




!returns the inflection point of the potential
  function ncli_x_inflection(alpha,phi0,n)
    implicit none
    real(kp) :: ncli_x_inflection
    real(kp), intent(in) :: alpha,phi0,n

    ncli_x_inflection = phi0*(alpha/((3._kp+2._kp*n)*(4._kp+2._kp*n)))**(1._kp/(4._kp+2._kp*n))

  end function  ncli_x_inflection



!return the value of mu at which the inflection point and the
!vanishing of the potential coincide. For values of phi0 smaller than
!phizeromin, the inflection point and the plateau would reside in a region in
!which the potential is negative which is not acceptable
  function ncli_phizeromin(alpha,n)
    implicit none
    real(kp) :: ncli_phizeromin
    real(kp), intent(in) :: alpha, n

    ncli_phizeromin = (sqrt(alpha/(4._kp+2._kp*n))*exp((2._kp + n)/alpha) &
         *sqrt(exp(1._kp/(3._kp+2._kp*n))/(3._kp+2._kp*n)))**(-1._kp/(2._kp+n))

  end function ncli_phizeromin


! returns the location of the two zeros of eps2
  function ncli_x_epstwozero(alpha,phi0,n)
    implicit none
    real(kp), dimension(2) :: ncli_x_epstwozero
    real(kp), intent(in) :: alpha, phi0,n
    real(kp), parameter :: tolFind=tolkp

    real(kp), dimension(2), save :: xeps2 = -1._kp
    real(kp), save :: alphasav=-1._kp, phi0sav=-1._kp, nsav=-1._kp
!$omp threadprivate(alphasav,phi0sav,nsav,xeps2)

    type(transfert) :: ncliData
    real(kp) :: xguess, eps2, newton, mini, maxi
    real(kp) :: x, V, dV, eps3

    integer :: count
    integer, parameter :: countMax = 10._kp

    if (phi0.lt.ncli_phizeromin(alpha,n)) stop 'ncli_x_epsoneunity: phi0 < phi0min!'

    if ((alpha.eq.alphasav).and.(phi0.eq.phi0sav).and.(n.eq.nsav)) then
       ncli_x_epstwozero = xeps2
       return
    else
       alphasav = alpha
       phi0sav = phi0
       nsav = n
    endif


    xguess = ncli_x_inflection(alpha,phi0,n)

    eps2 = 1._kp
    x = xguess
    count = 0

!Newton's method to find root closest to the inflection point, which
!is the smallest root of eps=2
    do while (abs(eps2).gt.tolFind)

!eps2/eps2'
       newton = (2*(x**(5 + 2*n) + phi0**(4 + 2*n)*x*(1 + alpha*Log(x)))*(2*(2 + n)*x**(8 + 4*n) &
            + alpha*phi0**(8 + 4*n)*(1 + alpha + alpha*Log(x)) &
            + phi0**(4 + 2*n)*x**(4 + 2*n)*(-2*(2 + n)*(3 + 2*n) + alpha*(9 + 4*n) &
            - 2*alpha*(2 + n)*(3 + 2*n)*Log(x))))/((x**(5 + 2*n) + phi0**(4 + 2*n)*x*(1 &
            + alpha*Log(x)))*((2*alpha**2*phi0**(8 + 4*n))/x + 16*(2 + n)**2*x**(7 + 4*n) &
            - 8*(2 + n)*phi0**(4 + 2*n)*x**(3 + 2*n)*(-(alpha*(3 + n)) + (2 + n)*(3 + 2*n) &
            + alpha*(2 + n)*(3 + 2*n)*Log(x))) - 2*((5 + 2*n)*x**(4 + 2*n) + phi0**(4 + 2*n)*(1 &
            + alpha + alpha*Log(x)))*(4*(2 + n)*x**(8 + 4*n) + 2*alpha*phi0**(8 + 4*n)*(1 &
            + alpha + alpha*Log(x)) - 2*phi0**(4 + 2*n)*x**(4 + 2*n)*(2*(2 + n)*(3 + 2*n) &
            - alpha*(9 + 4*n) + 2*alpha*(2 + n)*(3 + 2*n)*Log(x))))

!which is also (less robust)
!       eps3 = ncli_epsilon_three(x,alpha,phi0,n)
!       V = ncli_norm_potential(x,alpha,phi0,n)
!       dV = ncli_norm_deriv_potential(x,alpha,phi0,n)
!       newton = - dV/V/eps3


       x = x - newton


       eps2 = ncli_epsilon_two(x,alpha,phi0,n)

!numerical limitation
       if (abs(newton).lt.epsilon(1._kp)) exit

!safeguard
       count = count + 1
       if (count.eq.countMax) then
          write(*,*)'alpha= phi0= n= ',alpha,phi0,n
          write(*,*)'x= eps2= ',x,eps2
          write(*,*)'Iterated= ',countMax
          stop 'ncli_x_epstwozero: failed to converge'
       endif
    enddo

    xeps2(1) = x

!Then zbrenting the other one

    ncliData%real1 = alpha
    ncliData%real2 = phi0
    ncliData%real3 = n

    mini = x * (1._kp + 10._kp*epsilon(1._kp))
    maxi = xNumAccMax

    xeps2(2) = zbrent(find_ncli_x_epstwozero,mini,maxi,tolFind,ncliData)

    ncli_x_epstwozero = xeps2

  end function ncli_x_epstwozero



  function find_ncli_x_epstwozero(x,ncliData)
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ncliData
    real(kp) :: find_ncli_x_epstwozero
    real(kp) :: alpha,phi0,n

    alpha = ncliData%real1
    phi0 = ncliData%real2
    n = ncliData%real3

    find_ncli_x_epstwozero = ncli_epsilon_two(x,alpha,phi0,n)

  end function find_ncli_x_epstwozero



!returns the value of eps1 in its first minimum
  function ncli_epsilon_one_min(alpha,phi0,n)
    implicit none
    real(kp) :: ncli_epsilon_one_min
    real(kp), intent(in) :: alpha, phi0, n
    real(kp), dimension(2) :: xeps2

    xeps2 = ncli_x_epstwozero(alpha,phi0,n)
    ncli_epsilon_one_min = ncli_epsilon_one(xeps2(1),alpha,phi0,n)

  end function ncli_epsilon_one_min


!returns the value of eps1 at its first maximum
  function ncli_epsilon_one_max(alpha,phi0,n)
    implicit none
    real(kp) :: ncli_epsilon_one_max
    real(kp), intent(in) :: alpha, phi0, n
    real(kp), dimension(2) :: xeps2

    xeps2 = ncli_x_epstwozero(alpha,phi0,n)
    ncli_epsilon_one_max = ncli_epsilon_one(xeps2(2),alpha,phi0,n)

  end function ncli_epsilon_one_max



!returns a vector containing the three roots of eps1=1
  function ncli_x_epsoneunity(alpha,phi0,n)
    implicit none
    real(kp), dimension(3) :: ncli_x_epsoneunity
    real(kp), intent(in) :: alpha, phi0, n
    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: ncliData

    real(kp), dimension(3), save :: xeps1 = -1._kp
    real(kp), save :: alphasav=-1._kp, phi0sav=-1._kp, nsav=-1._kp
!$omp threadprivate(alphasav,phi0sav,nsav,xeps1)

    logical :: oneroot
    real(kp), dimension(2) ::  xeps2
    real(kp) :: eps1min, eps1max, xzero
    real(kp) :: mini, maxi

    if (phi0.lt.ncli_phizeromin(alpha,n)) stop 'ncli_x_epsoneunity: phi0 < phi0min!'

    if ((alpha.eq.alphasav).and.(phi0.eq.phi0sav).and.(n.eq.nsav)) then
       ncli_x_epsoneunity = xeps1
       return
    else
       alphasav = alpha
       phi0sav = phi0
       nsav = n
    endif

    xzero = ncli_x_potzero(alpha,phi0,n)
    xeps2 = ncli_x_epstwozero(alpha,phi0,n)
    eps1min = ncli_epsilon_one(xeps2(1),alpha,phi0,n)
    eps1max = ncli_epsilon_one(xeps2(2),alpha,phi0,n)

    ncliData%real1 = alpha
    ncliData%real2 = phi0
    ncliData%real3 = n

!eps1->infty in 0, there is only one root in these cases only
    oneroot = ((eps1min.gt.1._kp).and.(eps1max.gt.1._kp)) &
         .or. ((eps1min.lt.1._kp).and.(eps1max.lt.1._kp))

    if (oneroot) then

       mini = xzero + epsilon(1._kp)
       maxi = xNumAccMax

       xeps1(:) = zbrent(find_ncli_x_epsoneunity,mini,maxi,tolFind,ncliData)

!tworoots
    elseif ((eps1min.eq.1._kp).and.(eps1max.gt.1._kp)) then

       xeps1(1) = xeps2(1)
       xeps1(2) = xeps2(1)

       mini = xeps2(2)
       maxi = xNumAccMax

       xeps1(3) = zbrent(find_ncli_x_epsoneunity,mini,maxi,tolFind,ncliData)

!tworoots
    elseif ((eps1min.lt.1._kp).and.(eps1max.eq.1._kp)) then

       xeps1(3) = xeps2(2)

       mini = xzero + epsilon(1._kp)
       maxi = xeps2(2)

       xeps1(1) = zbrent(find_ncli_x_epsoneunity,mini,maxi,tolFind,ncliData)
       xeps1(2) = xeps1(1)

!threeroots
    else
       mini = xzero + epsilon(1._kp)
       maxi = xeps2(1)
       xeps1(1) = zbrent(find_ncli_x_epsoneunity,mini,maxi,tolFind,ncliData)

       mini = xeps2(1)
       maxi = xeps2(2)
       xeps1(2) = zbrent(find_ncli_x_epsoneunity,mini,maxi,tolFind,ncliData)

       mini = xeps2(2)
       maxi = xNumAccMax
       xeps1(3) = zbrent(find_ncli_x_epsoneunity,mini,maxi,tolFind,ncliData)

    endif

    ncli_x_epsoneunity(1:3) = xeps1(1:3)


  end function ncli_x_epsoneunity



  function find_ncli_x_epsoneunity(x,ncliData)
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ncliData
    real(kp) :: find_ncli_x_epsoneunity
    real(kp) :: alpha,phi0,n

    alpha = ncliData%real1
    phi0 = ncliData%real2
    n = ncliData%real3

    find_ncli_x_epsoneunity = alpha*(phi0*phi0)**(2+n) &
         + 2._kp*(2._kp+n)*(x*x)**(2._kp+n) - sqrt(2._kp)*x &
         * ((x*x)**(2+n) + (phi0*phi0)**(2+n)*(1._kp+alpha*log(x)))


  end function find_ncli_x_epsoneunity




!returns x at the end of inflation defined as epsilon1=1
  function ncli_x_endinf(alpha,phi0,n)
    implicit none
    real(kp), intent(in) :: alpha,phi0,n
    real(kp) :: ncli_x_endinf
    real(kp), parameter :: tolFind=tolkp

    real(kp), dimension(3) :: xeps1

    xeps1 = ncli_x_epsoneunity(alpha,phi0,n)

!the smallest value
    ncli_x_endinf = xeps1(1)

  end function ncli_x_endinf




!returns the maximum value for xini such that slow-roll inflation is
!valid (the trajectory should not cross the xeps1(3) value if this one
!exists)
  function ncli_xinimax(alpha,phi0,n)
    implicit none
    real(kp), intent(in) :: alpha,phi0,n
    real(kp) :: ncli_xinimax
    real(kp), parameter :: tolFind=tolkp
    real(kp), dimension(3) :: xeps1

    xeps1 = ncli_x_epsoneunity(alpha,phi0,n)

!happens only when eps1max > 1
    if (xeps1(2).ne.xeps1(3)) then
       ncli_xinimax = xeps1(2)
    else
       ncli_xinimax = xNumAccMAx
    endif

  end function ncli_xinimax



!returns the minimal value of xini such that there is at least efold
!slow-roll inflation
  function ncli_xinimin(efold,alpha,phi0,n)
    implicit none
    real(kp) :: ncli_xinimin
    real(kp), intent(in) :: efold,alpha,phi0,n
    real(kp) :: xiniMax, xend, efoldMax
    logical, parameter :: display = .false.

    xiniMax = ncli_xinimax(alpha,phi0,n)
    xend = ncli_x_endinf(alpha,phi0,n)

    efoldMax = -ncli_efold_primitive(xend,alpha,phi0,n) &
         +ncli_efold_primitive(xiniMax,alpha,phi0,n)


    if (efold.gt.efoldMax) then
       if (display) then
          write(*,*)'ncli_xinimin: not enough efolds!'
          write(*,*)'alpha= phi0= n= ',alpha,phi0,n
          write(*,*)'efold requested= ',efold,'efold maxi= ',efoldMax
       endif
       ncli_xinimin = xiniMax
       return
    endif

    ncli_xinimin = ncli_x_trajectory(efold,xiniMax,alpha,phi0,n)

  end function ncli_xinimin



  !this is integral[V(phi)/V'(phi) dphi]
  function ncli_efold_primitive(x,alpha,phi0,n)
    implicit none
    real(kp), intent(in) :: x,alpha,phi0,n
    real(kp) :: ncli_efold_primitive
    type(transfert) :: ncliData
!avoids prohibitive integration time in QUADPREC
    real(kp), parameter :: tolInt = max(tolkp,toldp)

    integer, parameter :: neq = 1
    real(kp) :: xvar, xinf
    real(kp), dimension(neq) :: yvar

!let us start at the inflection point
    xvar = ncli_x_inflection(alpha,phi0,n)
    yvar(1) = 0._kp

    ncliData%real1 = alpha
    ncliData%real2 = phi0
    ncliData%real3 = n

    call easydverk(neq,find_ncli_efold_primitive,xvar,yvar,x,tolInt,ncliData)

    ncli_efold_primitive = yvar(1)

  end function ncli_efold_primitive

  subroutine find_ncli_efold_primitive(n,x,y,yprime,ncliData)
    implicit none
    integer :: n
    real(kp) :: x
    real(kp), dimension(n) :: y, yprime
    type(transfert), optional, intent(inout) :: ncliData
    real(kp) :: alpha, phi0, nncli, x4p2n, phi4p2n

    alpha = ncliData%real1
    phi0 = ncliData%real2
    nncli = ncliData%real3

    x4p2n = (x*x)**(2._kp+nncli)
    phi4p2n = (phi0*phi0)**(2._kp+nncli)

    yprime(1) = x*(x4p2n + phi4p2n*(1._kp+alpha*log(x))) &
         / (alpha*phi4p2n + 2._kp*(2._kp+nncli)*x4p2n)

!    print *,'yprime',yprime(1)
!regularized to avoid eps1=0
!    print *,'test', 1._kp/sqrt(epsilon(1._kp) + 2._kp*ncli_epsilon_one(x,alpha,phi0,nncli))


  end subroutine find_ncli_efold_primitive


  !returns x at bfold=-efolds before the end of inflation
  function ncli_x_trajectory(bfold,xend,alpha,phi0,n)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,phi0,n
    real(kp) :: ncli_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ncliData

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = xNumAccMax

    ncliData%real1 = alpha
    ncliData%real2 = phi0
    ncliData%real3 = n
    ncliData%real4 = -bfold + ncli_efold_primitive(xend,alpha,phi0,n)

    ncli_x_trajectory = zbrent(find_ncli_x_trajectory,mini,maxi,tolFind,ncliData)

  end function ncli_x_trajectory

  function find_ncli_x_trajectory(x,ncliData)
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ncliData
    real(kp) :: find_ncli_x_trajectory
    real(kp) :: alpha,phi0,n,NplusPrimEnd

    alpha=ncliData%real1
    phi0=ncliData%real2
    n=ncliData%real3
    NplusPrimEnd = ncliData%real4

    find_ncli_x_trajectory = ncli_efold_primitive(x,alpha,phi0,n) - NplusPrimEnd

  end function find_ncli_x_trajectory



end module nclisr
